
# ----------------------------------------------------------------------------- #
# R script for the random forest model from the journal article
# "Quantifying mercury stocks and fluxes from permafrost coastal erosion using machine learning"
# by Jaspers, K., Irrgang, A.M., Wolter, J., Haugk, C., Petzold, P.,
# Jonsson, S., Lantuit, H., Fritz, M.
# published in Communications Earth & Environment
#
# Script name: random_forest_model.R
#
# Purpose of this script: This script provides the code to the random forest
# model and its outcomes that are published in the above-mentioned article.
#
# Author: Katharina Jaspers
# Last updated: July 21, 2026
# Copyright (c) Katharina Jaspers, 2026
# Contact: katharina.jaspers@awi.de
#
# If you use this code, please cite:
# Jaspers, K., Irrgang, A.M., Wolter, J., Haugk, C., Petzold, P., Jonsson, S., Lantuit, H., Fritz, M. (2026): 
# Stocks, fluxes, and fate: quantifying mercury release from permafrost coastal erosion. 
# Communications Earth & Environment.
#
# ----------------------------------------------------------------------------- #
#
# To run this script, you need the following files in the same folder as this
# script, or in the folder defined below as project_dir:
# - all_layers_data.xlsx
# - terrain_units.gpkg
# - landcover_raster.tif
#
# List of abbreviations:
# - hg = mercury
# - imp = importance
# - lc = land cover
# - lm = linear model
# - mehg = methylmercury
# - miss = missing
# - MLR = multiple linear regression
# - mod = modeled/to be modeled
# - obs = observed
# - pred = predicted/to be predicted
# - rf = random forest
# - train = training
# - tu = terrain unit
# - valid = validation
# ----------------------------------------------------------------------------- #







# 0 setup ######################################################################

## values
random_seed <- 123
n_tuning_splits <- 30
n_model_runs <- 500
n_tuning_trees <- 500
n_rf_trees <- 500
boruta_max_runs <- 200
boruta_p_value <- 0.05

# Permutation importance is not used by any subsequent calculation. Set this to
# TRUE only if the (computationally expensive) per-run iml objects are needed.
calculate_feature_importance <- FALSE

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
set.seed(random_seed)

## load libraries
required_packages <- c(
  "Boruta",
  "caret",
  "dplyr",
  "forcats",
  "iml",
  "progressr",
  "quantregForest",
  "sf",
  "terra",
  "pangaear"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the following package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))


## paths
# The default assumes that all input files are stored next to this script.
get_script_dir <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = TRUE)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

# set working directory
project_dir <- get_script_dir()

# set paths to files
input_paths <- list(
  terrain_units = file.path(project_dir, "terrain_units.gpkg"),
  land_cover = file.path(project_dir, "landcover_raster.tif")
)

# check for missing files
missing_input <- names(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing_input) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_input, collapse = ", "),
    call. = FALSE
  )
}


## read and setup Hg data
# load data from PANGAEA
# Jaspers, Katharina; Couture, Nicole; Wolter, Juliane; Fritz, Michael (2026): 
# Mercury (Hg) contents for terrain units along the Yukon Coastal Plain, Canada, including geomorphological parameters [dataset]. 
# PANGAEA, https://doi.org/10.1594/PANGAEA.987813

file <- pg_data(doi = "10.1594/PANGAEA.987813")
hg_data <- file[[1]]$data

# delete columns that are not needed
hg_data <- hg_data[, -which(names(hg_data) %in% c("Event", "Device (Device/method that was used t...)", 
                                                  "Sampling date (Date of sampling of the Hg sa...)", 
                                                  "Campaign (Project during which the Hg s...)"))]

# rename columns
col_mapping <- c(
  "sample_no"    = "Sample no",
  "sample_id"    = "Sample ID (Sample ID from Couture et al....)",
  "tu"           = "Location (Terrain unit)",
  "latitude"     = "Latitude",
  "longitude"    = "Longitude",
  "alt"          = "ALD [cm] (Active layer thickness)",
  "al_p"         = "Comment (Active layer or permafrost? [...)",
  "depth_min"    = "Depth soil min [m]",
  "depth_max"    = "Depth soil max [m]",
  "depth_mean"   = "Depth soil [m] (mean)",
  "lithology"    = "Lithology (Surficial lithology)",
  "TOC"          = "TOC [%]",
  "ref_TOC"      = "Reference (Source for TOC value)",
  "DBD"          = "DBD [g/cm**3]",
  "ref_DBD"      = "Reference (Source for DBD value)",
  "Hg"           = "Hg [µg/kg]",
  "Hg_sd"        = "Hg std dev [±] (calculated from repeated meas...)",
  "Hg_ref"       = "Reference (Source for Hg value)",
  "grain_size"   = "Grain size descr (Qualitative grain size descri...)",
  "slope_score"  = "Slope (Main slope score of the respe...)",
  "bluff_height" = "Bluff h [m]",
  "change_rate"  = "Change rate [m/a]",
  "area_eroded"  = "Area erod [m**2/a]",
  "ice_volume"   = "Ice [%] (Massive and wedge ice)"
)

# check for missing columns
missing_in_data <- setdiff(col_mapping, colnames(hg_data))
missing_in_mapping <- setdiff(colnames(hg_data), col_mapping)

if (length(missing_in_data) > 0) {
  stop("Expected columns missing in the dataset: ", paste(missing_in_data, collapse = ", "))
}
if (length(missing_in_mapping) > 0) {
  warning("Columns in the dataset that are not considered in the code: ", paste(missing_in_mapping, collapse = ", "))
}

hg_data <- dplyr::rename(hg_data, all_of(col_mapping))

# convert column types
hg_data <- hg_data %>%
  dplyr::mutate(across(where(is.character), as.factor),
                slope_score = as.factor(slope_score))


## read and setup spatial layer of terrain units
tu_sf <- st_read(input_paths$terrain_units)
tu_sf <- tu_sf %>%
  dplyr::select(any_of(c("Segment_na", "GEOLUNIT", "Unit_No", "area_terr", "Lc", "geom")))
names(tu_sf) <- c("tu", "geol_unit", "unit_no", "area_terr", "Lc", "geom")

# sort by unit_no
tu_sf <- tu_sf[order(tu_sf$unit_no, decreasing = TRUE), ]


## read and setup land cover raster
lc_rast <- rast(input_paths$land_cover)

## calculate the most common land cover classes per terrain unit
# count number of cells per land cover class and tu
lc_zonal <- terra::extract(lc_rast, tu_sf)
lc_zonal <- na.omit(lc_zonal)

# convert IDs (created by terra::extract) to tu
id_to_tu <- data.frame(ID = 1:nrow(tu_sf), tu = tu_sf$tu)
lc_zonal <- left_join(lc_zonal, id_to_tu, by = "ID")

# split raster values into different tables per tu
lc_by_tu <- split(lc_zonal$landcover_raster, lc_zonal$tu)

# function to get top five lc classes and their spatial proportion per terrain unit
get_top5_lc <- function(class_values) {
  tab <- table(class_values)
  total <- sum(tab)
  sorted <- sort(tab, decreasing = TRUE)
  
  # get top five
  top_5 <- head(sorted, 5)
  classes <- as.integer(names(top_5))
  percents <- round(100 * as.numeric(top_5) / total, 1)
  
  top5_names <- c()
  for (i in 1:5) {
    top5_names <- c(top5_names, paste0("top", i, "_class"), paste0("top", i, "_percent"))
  }
  
  result <- as.list(c(rbind(classes, percents)))
  names(result) <- top5_names
  return(result)
}

# apply function
lc_top5_by_tu <- lapply(lc_by_tu, get_top5_lc)

lc_top5 <- bind_rows(lc_top5_by_tu, .id = "tu")

lc_top5 <- merge(tu_sf, lc_top5, by.x = "tu", by.y = "tu", all.x = TRUE)

# attach land cover values to hg_data
hg_data <- hg_data %>%
  left_join(
    lc_top5 %>%
      st_drop_geometry() %>%
      dplyr::select(
        tu, top1_class, top1_percent, top2_class, top2_percent, top3_class, 
        top3_percent, top4_class, top4_percent, top5_class, top5_percent
      ),
    by = "tu"
  )

# -> hg_data now has ten extra columns with the top 5 land cover classes 
# and their spatial coverage in %





# 1 data preparation ####################################

# create data set with (potential) predictors for random forest model
rf_data <- hg_data %>% dplyr::select(c(sample_no, tu, latitude, longitude, al_p, depth_mean, lithology, 
                                       TOC, DBD, Hg, Hg_sd, grain_size, slope_score, 
                                       top1_class, top1_percent, top2_class, top2_percent,
                                       top3_class, top3_percent, top4_class, top4_percent, 
                                       top5_class, top5_percent, depth_min, depth_max)) %>% 
  dplyr::mutate(tu = as.factor(tu), 
                top1_class = as.factor(top1_class),
                top2_class = as.factor(top2_class),
                top3_class = as.factor(top3_class),
                top4_class = as.factor(top4_class),
                top5_class = as.factor(top5_class))

# create subset: only rows where Hg values are missing
rf_miss <- rf_data[which(is.na(rf_data$Hg)),]

# create subset: only rows where Hg values are available
rf_obs <- rf_data[which(!is.na(rf_data$Hg)),]





# 2 RF model & prediction ####################################

## define relevant functions

# create training-data IDs stratified by terrain unit 
create_train_ids <- function(data) {
  data %>%
    dplyr::group_by(tu) %>%
    group_modify(~ {
      n_rows <- nrow(.x)
      n_sample <- floor(0.75 * n_rows)
      slice_sample(.x, n = n_sample)
    }) %>%
    ungroup() %>%
    pull(sample_no)
}

# match validation/prediction factor levels to the training data
# the rarest level in the training dataset is set to "Other"
# factor levels present in other but absent from train are set to "Other"
harmonize_factor_levels <- function(train, other, columns) {
  for (col in columns) {
    
    train_values <- unique(as.character(train[[col]]))
    other_values <- unique(as.character(other[[col]]))
    
    values_only_in_other <- setdiff(other_values, train_values)
    
    if (length(values_only_in_other) > 0) {
      
      train[[col]] <- fct_lump(
        train[[col]],
        n = length(unique(train[[col]])) - 1
      )
      
      other[[col]] <- fct_other(
        other[[col]],
        keep = levels(train[[col]])
      )
      
      stopifnot(
        identical(
          levels(train[[col]]),
          levels(other[[col]])
        )
      )
    }
  }
  
  list(train = train, other = other)
}

# function to calculate Root Mean Squared Error (RMSE)
calc_rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

# function to calculate Mean Absolute Error (MAE)
calc_mae <- function(obs, pred) {
  mean(abs(obs - pred))
}

# function to calculate coefficient of determination (R^2)
calc_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) /
    sum((obs - mean(obs))^2)
}



## 2.1 tuning ####

# tune a small parameter grid over several random train-validation splits
# the aim is not exhaustive optimization, but a robust parameter choice

set.seed(random_seed)

rf_tune_grid <- expand.grid(
  mtry = c(2,4,6,8,10,12,14,16,18),
  nodesize = c(3, 5, 7),
  ntree = n_tuning_trees
)

rf_tune_runs <- vector("list", n_tuning_splits * nrow(rf_tune_grid))
counter <- 1

for(i in 1:n_tuning_splits) {
  
  set.seed(i*i)
  
  # create training data set ids
  rf_train_ids <- create_train_ids(rf_obs)
  
  # create training and validation data sets
  rf_train_full <- rf_obs %>%
    dplyr::filter(sample_no %in% rf_train_ids)
  
  rf_valid_full <- rf_obs %>%
    dplyr::filter(!sample_no %in% rf_train_ids)
  
  # factor columns except for tu which will not serve as a predictor variable
  fac_cols <- names(rf_train_full)[-which(names(rf_train_full) == "tu")]
  fac_cols <- fac_cols[sapply(rf_train_full[fac_cols], is.factor)]
  
  # harmonize factor levels
  harmonized <- harmonize_factor_levels(rf_train_full, rf_valid_full, fac_cols)
  rf_train_full <- harmonized$train
  rf_valid_full <- harmonized$other
  
  # Boruta variable selection
  rf_boruta_df <- rf_train_full[, -which(names(rf_train_full) %in% c("sample_no", "tu", "Hg_sd"))]
  
  set.seed(i*i)
  
  boruta_output <- Boruta(
    Hg ~ .,
    data = rf_boruta_df,
    doTrace = 0,
    pValue = boruta_p_value,
    maxRuns = boruta_max_runs
  )
  
  selected_vars <- getSelectedAttributes(boruta_output, withTentative = TRUE)
  
  rf_train_tune <- rf_train_full[, c("Hg", selected_vars)]
  rf_valid_tune <- rf_valid_full[, c("sample_no", "Hg", selected_vars)]
  
  x_train <- rf_train_tune[, names(rf_train_tune) != "Hg"]
  y_train <- rf_train_tune$Hg
  
  x_valid <- rf_valid_tune[, !(names(rf_valid_tune) %in% c("sample_no", "Hg"))]
  y_valid <- rf_valid_tune$Hg
  
  for(j in seq_len(nrow(rf_tune_grid))) {
    
    params <- rf_tune_grid[j, ]
    
    # mtry must not be larger than the number of selected predictors
    current_mtry <- min(params$mtry, ncol(x_train))
    
    set.seed(j * 1000 + i)
    
    rf_tune_fit <- quantregForest(
      x = x_train,
      y = y_train,
      ntree = params$ntree,
      mtry = current_mtry,
      nodesize = params$nodesize,
      keep.inbag = TRUE
    )
    
    pred <- predict(
      rf_tune_fit,
      newdata = x_valid,
      what = c(0.025, 0.5, 0.975))
    
    if(is.null(pred)) next
    
    rf_pred_tune <- as.data.frame(pred)
    
    pred_low <- rf_pred_tune$`quantile= 0.025`
    pred_med <- rf_pred_tune$`quantile= 0.5`
    pred_upp <- rf_pred_tune$`quantile= 0.975`
    
    rf_tune_runs[[counter]] <- data.frame(
      run = i,
      mtry = params$mtry,
      mtry_used = current_mtry,
      nodesize = params$nodesize,
      ntree = params$ntree,
      n_predictors = ncol(x_train),
      RMSE = sqrt(mean((y_valid - pred_med)^2, na.rm = TRUE)),
      MAE = mean(abs(y_valid - pred_med), na.rm = TRUE),
      bias = mean(pred_med - y_valid, na.rm = TRUE),
      coverage_95 = mean(y_valid >= pred_low & y_valid <= pred_upp, na.rm = TRUE)
    )
    
    counter <- counter + 1
  }
  
  print(paste("Tuning split", i, "of", n_tuning_splits, "done"))
}

rf_tune_results <- bind_rows(rf_tune_runs)

rf_tune_summary <- rf_tune_results %>%
  dplyr::group_by(mtry, nodesize, ntree) %>%
  dplyr::summarise(
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_MAE = mean(MAE),
    sd_MAE = sd(MAE),
    mean_bias = mean(bias),
    mean_coverage_95 = mean(coverage_95),
    .groups = "drop"
  ) %>%
  arrange(mean_RMSE)

rf_tune_summary

rf_best_params <- rf_tune_summary[1, ]



## 2.2 training & evaluation ####
# quantile regression forest (qrf) model, following Vaysse & Lagacherie 2017
# train and predict the model repeatedly, each time varying the training data set
# repeat model training to capture variability in predictions

handlers("txtprogressbar")
handlers(global = FALSE)

## training and evaluation
rf_run_results <- with_progress({
  
  p <- progressor(steps = n_model_runs)
  
  lapply(
    seq_len(n_model_runs),
    function(i) {
      
      set.seed(i)
      
      # create training data set ids
      rf_train_ids <- create_train_ids(rf_obs)
      
      # create training data set
      rf_train <- rf_obs %>%
        dplyr::filter(sample_no %in% rf_train_ids)
      
      # create validation data set (all values that were not used for training)
      rf_val <- rf_obs %>%
        dplyr::filter(!sample_no %in% rf_train_ids)
      
      # harmonize factor levels
      fac_cols <- names(rf_train)[-which(names(rf_train) == "tu")]
      fac_cols <- fac_cols[sapply(rf_train[fac_cols], is.factor)]
      
      harmonized <- harmonize_factor_levels(rf_train, rf_val, fac_cols)
      rf_train <- harmonized$train
      rf_val <- harmonized$other
      
      # Boruta
      rf_boruta_df <- rf_train[, -which(names(rf_train) %in% c("sample_no", "tu", "Hg_sd"))]
      set.seed(i)
      boruta_output <- Boruta(Hg ~ ., data = rf_boruta_df, doTrace=0, pValue = boruta_p_value,
                              maxRuns = boruta_max_runs)
      
      selected_features <- getSelectedAttributes(
        boruta_output,
        withTentative = TRUE
      )
      
      # keep only columns that were considered important by the Boruta algorithm
      rf_train <- rf_train[, c("Hg", selected_features)]
      rf_val   <- rf_val[, c("sample_no", "Hg", selected_features)]
      
      # create data frame for storing prediction results
      rf_run_pred <- data.frame("sample_no" = rf_val$sample_no,
                                "model_run" = i, 
                                "predicted" = NA, 
                                "Hg_conc_low" = NA, 
                                "Hg_conc_upp" = NA, 
                                "Hg" = NA)
      
      # predict Hg concentrations of the validation data and store them in a data frame
      set.seed(i * 1000)
      rf_mod <- quantregForest(x = rf_train[, names(rf_train) != "Hg"], 
                               y = rf_train$Hg, 
                               ntree = n_rf_trees, 
                               mtry = min(rf_best_params$mtry, ncol(rf_train[, names(rf_train) != "Hg"])), 
                               nodesize = rf_best_params$nodesize,
                               keep.inbag = TRUE)
      rf_quantiles <- predict(rf_mod, 
                              newdata = rf_val[, ! names(rf_val) %in% c("sample_no", "Hg")], 
                              what = c(0.025, 0.5, 0.975))
      
      rf_quantiles <- as.data.frame(rf_quantiles)
      
      rf_run_pred$predicted <- rf_quantiles$`quantile= 0.5`
      rf_run_pred$Hg_conc_upp <- rf_quantiles$`quantile= 0.975`
      rf_run_pred$Hg_conc_low <- rf_quantiles$`quantile= 0.025`
      rf_run_pred$Hg <- rf_val$Hg
      
      # quality check per run
      metrics <- data.frame(
        model_run = i,
        n_train = nrow(rf_train),
        n_valid = nrow(rf_val),
        
        RMSE = calc_rmse(rf_run_pred$Hg, rf_run_pred$predicted),
        MAE = calc_mae(rf_run_pred$Hg, rf_run_pred$predicted),
        R2 = calc_r2(rf_run_pred$Hg, rf_run_pred$predicted),
        r = cor(as.numeric(rf_run_pred$Hg),
                rf_run_pred$predicted, use = "complete.obs")
        
      )
      
      # calculate feature importance if set to TRUE
      importance <- NULL
      if (calculate_feature_importance) {
        predictor <- Predictor$new(
          model = rf_mod,
          data = rf_train[, -which(names(rf_train) == "Hg")],
          y = rf_train$Hg,
          predict.function = function(model, newdata) {
            predict(model, newdata = newdata, what = 0.5)
          }
        )
        
        set.seed(random_seed + i)
        importance <- FeatureImp$new(predictor, loss = "rmse")
      }
      
      result <- list(
        predictions = rf_run_pred,
        metrics = metrics,
        importance = importance,
        selected_features = selected_features
      )
      
      p(sprintf("Run %d", i))
      
      result
    }
  )
})

rf_metrics <- do.call(
  rbind,
  lapply(rf_run_results, `[[`, "metrics")
)

rf_pred <- do.call(
  rbind,
  lapply(rf_run_results, `[[`, "predictions")
)

# summarize outcomes from all model runs to get one value per data point
rf_pred_all <- rf_pred

rf_pred_all <- rf_pred_all %>%
  dplyr::group_by(sample_no) %>%
  dplyr::summarise(
    mean = mean(predicted),
    median = median(predicted),
    sd = sd(predicted),
    n_validated = n(),
    Hg_mean_upp = mean(Hg_conc_upp),
    Hg_mean_low = mean(Hg_conc_low),
    Hg = first(Hg),
    .groups = "drop"
  )


## quality measures averaged from all runs
rf_metrics_summary <- rf_metrics %>%
  dplyr::summarise(
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_MAE = mean(MAE),
    sd_MAE = sd(MAE),
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_r = mean(r),
    sd_r = sd(r)
  )

round(rf_metrics_summary,2)


## quality check overall

# correlation of actual and predicted Hg concentration
cor.test(rf_pred_all$Hg, rf_pred_all$mean)

# RMSE
calc_rmse(rf_pred_all$Hg, rf_pred_all$mean)

# sd
sd(rf_obs$Hg)

# R2
calc_r2(rf_pred_all$Hg, rf_pred_all$mean)

# MAE
calc_mae(rf_pred_all$Hg, rf_pred_all$mean)




## 2.3 prediction ####

# train a final model using all measured Hg data and predict Hg for the unmeasured data points

rf_final_train <- rf_obs
rf_final_miss <- rf_miss

## harmonize factor levels
fac_cols <- names(rf_final_train)[-which(names(rf_final_train) == "tu")]
fac_cols <- fac_cols[sapply(rf_final_train[fac_cols], is.factor)]

harmonized <- harmonize_factor_levels(
  rf_final_train,
  rf_final_miss,
  fac_cols
)

rf_final_train <- harmonized$train
rf_final_miss <- harmonized$other

## determine relevant predictors
rf_boruta_df <- rf_final_train[, -which(names(rf_final_train) %in% c("sample_no", "tu", "Hg_sd"))]

set.seed(random_seed)
boruta_output <- Boruta(Hg ~ ., 
                        data = rf_boruta_df, 
                        doTrace=2, 
                        pValue = boruta_p_value,
                        maxRuns = boruta_max_runs)

# see results from Boruta algorithm
getSelectedAttributes(boruta_output, withTentative = TRUE)
attStats(boruta_output)

# only keep relevant columns of training and prediction data sets
rf_final_train_sel <- rf_final_train[, c("Hg", 
                                         getSelectedAttributes(boruta_output, withTentative = TRUE))]
rf_final_miss_sel <- rf_final_miss[, c("Hg", 
                                       getSelectedAttributes(boruta_output, withTentative = TRUE))]

## one final model training and prediction
set.seed(random_seed)
rf_final_mod <- quantregForest(x = rf_final_train_sel[, -which(names(rf_final_train_sel) == "Hg")],
                               y = rf_final_train_sel$Hg, 
                               ntree = n_rf_trees, 
                               mtry = min(rf_best_params$mtry, ncol(rf_final_train_sel[, -which(names(rf_final_train_sel) == "Hg")])), 
                               nodesize = rf_best_params$nodesize, 
                               keep.inbag = TRUE)

rf_final_pred <- as.data.frame(predict(rf_final_mod, 
                                       newdata = rf_final_miss_sel[, -which(names(rf_final_miss_sel)=="Hg")],
                                       what = c(0.025, 0.5, 0.975)))
rf_final_pred <- cbind(rf_final_miss, rf_final_pred)

# -> rf_final_pred contains the predicted Hg concentrations for data points 
# where no Hg measurements were available





# 3 MLR ####

## 3.1 predictor determination ####

set.seed(random_seed)

# use the same measured dataset as for RF
lm_data <- rf_obs

## predictors for LM
# remove variables that should not be predictors
lm_candidate_vars <- setdiff(
  names(lm_data),
  c("sample_no", "tu", "Hg", "Hg_sd")
)

lm_data$Hg_log <- base::log(lm_data$Hg)

## define lm formula for factor evaluation
test_model <- lm(
  as.formula(
    paste(
      "Hg_log ~",
      paste(lm_candidate_vars, collapse=" + ")
    )
  ),
  data = lm_data
)

lm_alias <- alias(test_model)

# -> remove land cover class predictors to avoid redundant factor structures
# and unstable MLR coefficient estimates
lm_candidate_vars <- setdiff(
  lm_candidate_vars,
  c("top1_class", "top2_class", "top3_class", "top4_class", "top5_class")
)



## 3.2 evaluation ####

# storage
lm_raw_pred_runs <- vector("list", n_model_runs)
lm_log_pred_runs <- vector("list", n_model_runs)
lm_metric_runs <- vector("list", n_model_runs)

## raw (Hg) and log-transformed (log(Hg)) MLR and evaluation
for (i in seq_len(n_model_runs)) {
  
  set.seed(i)
  
  # same split logic as QRF: stratified by terrain unit
  lm_train_ids <- create_train_ids(lm_data)
  
  # prepare training and validation data
  lm_train <- lm_data %>%
    dplyr::filter(sample_no %in% lm_train_ids) %>%
    dplyr::select(sample_no, Hg, Hg_log, all_of(lm_candidate_vars))
  
  lm_val <- lm_data %>%
    dplyr::filter(!sample_no %in% lm_train_ids) %>%
    dplyr::select(sample_no, Hg, Hg_log, all_of(lm_candidate_vars))
  
  # determine correlated numeric predictors using training data only
  run_candidate_vars <- lm_candidate_vars
  
  num_vars <- run_candidate_vars[
    sapply(lm_train[run_candidate_vars], is.numeric)
  ]
  
  if (length(num_vars) > 1) {
    cor_matrix <- cor(
      lm_train[num_vars],
      use = "pairwise.complete.obs"
    )
    
    remove_cor <- findCorrelation(
      cor_matrix,
      cutoff = 0.85
    )
    
    run_candidate_vars <- setdiff(
      run_candidate_vars,
      num_vars[remove_cor]
    )
  }
  
  ## harmonize factor levels
  fac_cols <- names(lm_train)
  fac_cols <- fac_cols[sapply(lm_train[fac_cols], is.factor)]
  
  harmonized <- harmonize_factor_levels(lm_train, lm_val, fac_cols)
  lm_train <- harmonized$train
  lm_val <- harmonized$other
  
  # remove predictors with no variation in this training split
  usable_vars <- run_candidate_vars[
    sapply(lm_train[run_candidate_vars], function(x)
      length(unique(na.omit(x))) > 1
    )
  ]
  
  # MLR formula on raw Hg
  lm_formula_raw <- reformulate(
    termlabels = usable_vars,
    response = "Hg"
  )
  
  lm_formula_log <- reformulate(
    termlabels = usable_vars,
    response = "Hg_log"
  )
  
  # model building
  lm_raw_mod <- lm(lm_formula_raw, data = lm_train)
  lm_log_mod <- lm(lm_formula_log, data = lm_train)
  
  if(any(is.na(coef(lm_raw_mod)))) next
  if(any(is.na(coef(lm_log_mod)))) next
  
  ## predict validation data
  lm_pred_raw <- as.data.frame(
    predict(
      lm_raw_mod,
      newdata = lm_val,
      interval = "prediction",
      level = 0.95
    )
  )
  
  lm_pred_log <- as.data.frame(
    predict(
      lm_log_mod,
      newdata = lm_val,
      interval = "prediction",
      level = 0.95
    )
  )
  
  lm_raw_run_pred <- data.frame(
    sample_no = lm_val$sample_no,
    model_run = i,
    Hg = lm_val$Hg,
    predicted_raw = lm_pred_raw$fit,
    Hg_conc_low_raw = lm_pred_raw$lwr,
    Hg_conc_upp_raw = lm_pred_raw$upr
  )
  
  ## results from log-Hg model, both on log scale and back-transformed
  lm_log_run_pred <- data.frame(
    sample_no = lm_val$sample_no,
    model_run = i,
    Hg = lm_val$Hg,
    Hg_log = lm_val$Hg_log,
    predicted_log = lm_pred_log$fit,
    Hg_conc_low_log = lm_pred_log$lwr,
    Hg_conc_upp_log = lm_pred_log$upr,
    predicted_backtransformed = exp(lm_pred_log$fit),
    Hg_conc_low_backtransformed = exp(lm_pred_log$lwr),
    Hg_conc_upp_backtransformed = exp(lm_pred_log$upr)
  )
  
  lm_raw_pred_runs[[i]] <- lm_raw_run_pred
  lm_log_pred_runs[[i]] <- lm_log_run_pred
  
  ## fold-wise performance metrics
  lm_metric_runs[[i]] <- data.frame(
    model_run = i,
    n_train = nrow(lm_train),
    n_valid = nrow(lm_val),
    
    RMSE_raw = calc_rmse(lm_raw_run_pred$Hg, lm_raw_run_pred$predicted_raw),
    MAE_raw = calc_mae(lm_raw_run_pred$Hg, lm_raw_run_pred$predicted_raw),
    R2_raw = calc_r2(lm_raw_run_pred$Hg, lm_raw_run_pred$predicted_raw),
    r_raw = suppressWarnings(cor(lm_raw_run_pred$Hg,
                                 lm_raw_run_pred$predicted_raw,
                                 use = "complete.obs")),
    
    RMSE_log = calc_rmse(lm_log_run_pred$Hg_log, lm_log_run_pred$predicted_log),
    MAE_log = calc_mae(lm_log_run_pred$Hg_log, lm_log_run_pred$predicted_log),
    R2_log = calc_r2(lm_log_run_pred$Hg_log, lm_log_run_pred$predicted_log),
    r_log = suppressWarnings(cor(lm_log_run_pred$Hg_log,
                                 lm_log_run_pred$predicted_log,
                                 use = "complete.obs")),
    
    RMSE_backtransformed = calc_rmse(lm_log_run_pred$Hg, lm_log_run_pred$predicted_backtransformed),
    MAE_backtransformed = calc_mae(lm_log_run_pred$Hg, lm_log_run_pred$predicted_backtransformed),
    R2_backtransformed = calc_r2(lm_log_run_pred$Hg, lm_log_run_pred$predicted_backtransformed),
    r_backtransformed = suppressWarnings(cor(lm_log_run_pred$Hg,
                                             lm_log_run_pred$predicted_backtransformed,
                                             use = "complete.obs"))
  )
  
  if (i %% 50 == 0) {
    print(paste("LM run", i, "of", n_model_runs, "done"))
  }
  
}


# combine predictions and metrics
lm_raw_pred <- bind_rows(lm_raw_pred_runs)
lm_log_pred <- bind_rows(lm_log_pred_runs)
lm_metrics <- bind_rows(lm_metric_runs)

## summarize outcomes from all model runs to get one value per data point
lm_pred_all_raw <- lm_raw_pred %>%
  dplyr::group_by(sample_no) %>%
  dplyr::summarise(
    mean = mean(predicted_raw, na.rm = TRUE),
    median = median(predicted_raw, na.rm = TRUE),
    sd = sd(predicted_raw, na.rm = TRUE),
    n_validated = n(),
    Hg = first(Hg),
    .groups = "drop"
  )

lm_pred_all_log <- lm_log_pred %>%
  dplyr::group_by(sample_no) %>%
  dplyr::summarise(
    mean_log = mean(predicted_log, na.rm = TRUE),
    median_log = median(predicted_log, na.rm = TRUE),
    sd_log = sd(predicted_log, na.rm = TRUE),
    mean_backtransformed = mean(predicted_backtransformed, na.rm = TRUE),
    median_backtransformed = median(predicted_backtransformed, na.rm = TRUE),
    sd_backtransformed = sd(predicted_backtransformed, na.rm = TRUE),
    n_validated = n(),
    Hg = first(Hg),
    Hg_log = first(Hg_log),
    .groups = "drop"
  )

## quality measures averaged from all runs
lm_metrics_summary <- lm_metrics %>%
  dplyr::summarise(
    mean_RMSE_raw = mean(RMSE_raw, na.rm = TRUE),
    sd_RMSE_raw = sd(RMSE_raw, na.rm = TRUE),
    mean_MAE_raw = mean(MAE_raw, na.rm = TRUE),
    sd_MAE_raw = sd(MAE_raw, na.rm = TRUE),
    mean_R2_raw = mean(R2_raw, na.rm = TRUE),
    sd_R2_raw = sd(R2_raw, na.rm = TRUE),
    mean_r_raw = mean(r_raw, na.rm = TRUE),
    sd_r_raw = sd(r_raw, na.rm = TRUE),
    
    mean_RMSE_log = mean(RMSE_log, na.rm = TRUE),
    sd_RMSE_log = sd(RMSE_log, na.rm = TRUE),
    mean_MAE_log = mean(MAE_log, na.rm = TRUE),
    sd_MAE_log = sd(MAE_log, na.rm = TRUE),
    mean_R2_log = mean(R2_log, na.rm = TRUE),
    sd_R2_log = sd(R2_log, na.rm = TRUE),
    mean_r_log = mean(r_log, na.rm = TRUE),
    sd_r_log = sd(r_log, na.rm = TRUE),
    
    mean_RMSE_backtransformed = mean(RMSE_backtransformed, na.rm = TRUE),
    sd_RMSE_backtransformed = sd(RMSE_backtransformed, na.rm = TRUE),
    mean_MAE_backtransformed = mean(MAE_backtransformed, na.rm = TRUE),
    sd_MAE_backtransformed = sd(MAE_backtransformed, na.rm = TRUE),
    mean_R2_backtransformed = mean(R2_backtransformed, na.rm = TRUE),
    sd_R2_backtransformed = sd(R2_backtransformed, na.rm = TRUE),
    mean_r_backtransformed = mean(r_backtransformed, na.rm = TRUE),
    sd_r_backtransformed = sd(r_backtransformed, na.rm = TRUE)
  )

round(lm_metrics_summary, 2)


## quality check overall on raw Hg data

# correlation of actual and predicted Hg concentration
cor.test(lm_pred_all_raw$Hg, lm_pred_all_raw$mean)

# RMSE
calc_rmse(lm_pred_all_raw$Hg, lm_pred_all_raw$mean)

# R2
calc_r2(lm_pred_all_raw$Hg, lm_pred_all_raw$mean)

# MAE
calc_mae(lm_pred_all_raw$Hg, lm_pred_all_raw$mean)



## quality check overall on back-transformed log(Hg) data
# correlation
cor.test(lm_pred_all_log$Hg, lm_pred_all_log$mean_backtransformed)

# RMSE
calc_rmse(lm_pred_all_log$Hg, lm_pred_all_log$mean_backtransformed)

# R2
calc_r2(lm_pred_all_log$Hg, lm_pred_all_log$mean_backtransformed)

# MAE
calc_mae(lm_pred_all_log$Hg, lm_pred_all_log$mean_backtransformed)





# 4 stock and flux calculation ####

pred_layers <- rf_miss
obs_layers <- rf_obs

# insert Hg prediction intervals where Hg values had not been available
pred_layers$Hg <- rf_final_pred$`quantile= 0.5`
pred_layers$Hg_conc_low <- rf_final_pred$`quantile= 0.025`
pred_layers$Hg_conc_upp <- rf_final_pred$`quantile= 0.975`

# estimate prediction interval for measured Hg concentrations using the standard deviation
obs_layers$Hg_conc_low <- pmax(0, obs_layers$Hg - 2*obs_layers$Hg_sd)
obs_layers$Hg_conc_upp <- obs_layers$Hg + 2*obs_layers$Hg_sd

# combine predicted data with measured data to get a data set where all data points have a Hg value
all_layers <- bind_rows(pred_layers, obs_layers) %>%
  arrange(sample_no)

# add columns that were dropped for prediction before
all_layers <- all_layers %>%
  left_join(hg_data %>% dplyr::select(sample_no, sample_id, bluff_height, change_rate,
                                      area_eroded, ice_volume),
            by = "sample_no")

# convert tu to factor
all_layers$tu <- factor(all_layers$tu, levels = rev(unique(all_layers$tu)))

## derive MeHg concentrations
# the numbers represent the mean ratio of MeHg on Hg for active layer and permafrost, respectively
# derived from measured terrestrial samples (n = 52)
all_layers$mehg <- ifelse(all_layers$al_p == "al",
                          all_layers$Hg * 0.012,
                          all_layers$Hg * 0.003)

all_layers$mehg_conc_low <- ifelse(all_layers$al_p == "al",
                                   all_layers$Hg_conc_low * 0.012,
                                   all_layers$Hg_conc_low * 0.003)

all_layers$mehg_conc_upp <- ifelse(all_layers$al_p == "al",
                                   all_layers$Hg_conc_upp * 0.012,
                                   all_layers$Hg_conc_upp * 0.003)

## function to calculate Hg and MeHg stocks and fluxes
calculate_stocks_fluxes <- function(data,
                                    conc,
                                    conc_low = NULL,
                                    conc_upp = NULL,
                                    tu_df) {
  
  # concentration per layer
  data$stock_lay <- data$DBD * 1000 *
    (data$depth_max - data$depth_min) *
    (data[[conc]] * 1e-9)
  
  # correction for wedge ice
  data$stock_corr <- data$stock_lay -
    (data$stock_lay * data$ice_volume / 100)
  
  # terrain-unit summaries: contains one unique row per sample_no/layer, 
  # so layer stocks are summed within each terrain unit
  summary_df <- data %>%
    dplyr::group_by(tu) %>%
    dplyr::summarise(
      stock_m2 = sum(stock_corr),
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(
      tu_df %>%
        st_drop_geometry() %>%
        dplyr::select(tu, area_terr, Lc),
      by = "tu"
    ) %>%
    dplyr::mutate(
      stock_tu = stock_m2 * area_terr,
      flux_tu = stock_m2 * area_eroded * -sign(change_rate),
      flux_per_m = flux_tu / Lc
    )
  
  result <- list(
    data = data,
    summary = summary_df,
    total_stock = sum(summary_df$stock_tu),
    total_flux = sum(summary_df$flux_tu)
  )
  
  # uncertainty calculations
  if (!is.null(conc_low) && !is.null(conc_upp)) {
    
    data$stock_lay_low <- data$DBD * 1000 *
      (data$depth_max - data$depth_min) *
      (data[[conc_low]] * 1e-9)
    data$stock_corr_low <- data$stock_lay_low -
      (data$stock_lay_low * data$ice_volume / 100)
    
    data$stock_lay_upp <- data$DBD * 1000 *
      (data$depth_max - data$depth_min) *
      (data[[conc_upp]] * 1e-9)
    data$stock_corr_upp <- data$stock_lay_upp -
      (data$stock_lay_upp * data$ice_volume / 100)
    
    summary_uc <- data %>%
      dplyr::group_by(tu) %>%
      dplyr::summarise(
        stock_m2_low = sum(stock_corr_low),
        stock_m2_upp = sum(stock_corr_upp),
        across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      left_join(
        tu_df %>%
          st_drop_geometry() %>%
          dplyr::select(tu, area_terr, Lc),
        by = "tu"
      ) %>%
      dplyr::mutate(
        stock_tu_low = stock_m2_low * area_terr,
        stock_tu_upp = stock_m2_upp * area_terr,
        flux_tu_low = stock_m2_low * area_eroded * -sign(change_rate),
        flux_tu_upp = stock_m2_upp * area_eroded * -sign(change_rate)
      )
    
    result$uncertainty <- summary_uc
    result$total_stock_low <- sum(summary_uc$stock_tu_low)
    result$total_stock_upp <- sum(summary_uc$stock_tu_upp)
    result$total_flux_low <- sum(summary_uc$flux_tu_low)
    result$total_flux_upp <- sum(summary_uc$flux_tu_upp)
  }
  
  return(result)
}


## Hg results - complete YCP (surface elevation to sea floor)
hg_results <- calculate_stocks_fluxes(
  data = all_layers,
  conc = "Hg",
  conc_low = "Hg_conc_low",
  conc_upp = "Hg_conc_upp",
  tu_df = tu_sf
)

## upper 3 m of the YCP
upper3m_layers <- all_layers %>%
  ungroup() %>%
  dplyr::filter(depth_min < 3) %>%
  dplyr::mutate(depth_max = pmin(depth_max, 3))

upper3m_hg_results <- calculate_stocks_fluxes(
  data = upper3m_layers,
  conc = "Hg",
  conc_low = "Hg_conc_low",
  conc_upp = "Hg_conc_upp",
  tu_df = tu_sf
)

## print results
cat(paste(paste0("total Hg stock: ", round(hg_results$total_stock, -3), " kg (", 
                 round(hg_results$total_stock_low, -3), " to ", round(hg_results$total_stock_upp, -3), " kg)"), 
          paste0("total Hg flux: ", round(hg_results$total_flux, -1), " kg yr-1 (", 
                 round(hg_results$total_flux_low, -1), " to ", round(hg_results$total_flux_upp, -1), " kg yr-1)"), 
          paste0("upper 3 m Hg stock: ", round(upper3m_hg_results$total_stock, -3), " kg (", 
                 round(upper3m_hg_results$total_stock_low, -3), " to ", round(upper3m_hg_results$total_stock_upp, -3),
                 " kg)"), sep = "\n"), "\n"
)


## MeHg results - complete YCP (surface elevation to sea floor)
mehg_results <- calculate_stocks_fluxes(
  data = all_layers,
  conc = "mehg",
  conc_low = "mehg_conc_low",
  conc_upp = "mehg_conc_upp",
  tu_df = tu_sf
)

## print results
cat(paste(paste0("total MeHg stock: ", round(mehg_results$total_stock, -1), " kg (", 
                 round(mehg_results$total_stock_low, -1), " to ", round(mehg_results$total_stock_upp, -1), " kg)"), 
          paste0("total MeHg flux: ", round(mehg_results$total_flux, 2), " kg yr-1 (", 
                 round(mehg_results$total_flux_low, 2), " to ", round(mehg_results$total_flux_upp, 2), " kg yr-1)"), 
          sep = "\n"), "\n"
)



# End of script
