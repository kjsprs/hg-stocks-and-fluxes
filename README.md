# Quantifying mercury stocks and fluxes from permafrost coastal erosion using machine learning

This folder contains supplementary data to the journal article **Quantifying mercury stocks and fluxes from permafrost coastal erosion using machine learning**, published in Communications Earth & Environment.

## Repository Structure & Data Description
- `terrain_units.gpkg` - spatial layer of the study area and its terrain units, based on Couture et al. 2018 (https://doi.org/10.1002/2017JG004166)
- `LANDC_raster.tif` - spatial land cover raster of the study area, based on Bartsch et al. 2019 (https://doi.pangaea.de/10.1594/PANGAEA.897916)
- `random_forest_model.R` - main R script containing data processing, modeling, and prediction


## Requirements
The R script was created under R version 4.6.1

Required R packages are
- `Boruta`
- `caret`
- `dplyr`
- `forcats`
- `iml`
- `progressr`
- `quantregForest`
- `sf`
- `terra`
- `pangaear`


## How to Run
1. Download this repository
2. Open `random_forest_model.R` in R or RStudio
3. Set the working directory to the folder containing the repository files
4. Run the script in order. It will:
     - load required packages
     - read the data
     - process and prepare predictors
     - train and validate a random forest model
	 - train and validate a multiple linear regression model
     - predict Hg concentrations for unmeasured layers
     - calculate Hg stocks and fluxes
  

## Citation
If you use this code or data, please cite:

Jaspers, K., Irrgang, A.M., Wolter, J., Haugk, C., Petzold, P., Jonsson, S., Lantuit, H., Fritz. M. Quantifying mercury stocks and fluxes from permafrost coastal erosion using machine learning. Communications Earth & Environment (in press)


## Contact
For questions, please contact katharina.jaspers@awi.de or michael.fritz@awi.de

