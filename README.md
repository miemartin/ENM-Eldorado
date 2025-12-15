ENM-Eldorado

This repository contains the scripts used to build ecological niche models and estimate the potential seasonal distribution of Aedes aegypti and Aedes albopictus in Eldorado city, Misiones, Argentina.

The workflow integrates field occurrence data with remotely sensed environmental variables and applies presence-only ecological niche modeling using MaxEnt.

Software and environment

R (≥ 4.2)

QGIS 3.18.3 (used for initial remote sensing preprocessing)

Main R packages: terra, raster, sf, dismo, ENMeval, wallace, usdm, tidyverse

All spatial analyses are conducted in WGS84 for modeling and reprojected to UTM Zone 21S for area calculations.

Data overview
Occurrence data

Monthly ovitrap sampling at fixed sites from 2016 to 2018

Sites include tire shops, cemeteries, dwellings, and public green areas

Records were aggregated by season:

Spring–summer

Autumn–winter

A site was considered positive for a species in a season if at least one larva or pupa was detected

Repeated detections within the same site and season were not treated as independent records

Environmental predictors

Sentinel-2 derived indices: NDVI, NDWI, NDBI

Land cover classes: water, impervious surface, bare soil, herbaceous vegetation, tree cover

Distance-to rasters for water, impervious surface, and vegetation

Land Surface Temperature (LST) derived from Landsat 8

All predictors were resampled to 10 m resolution and cropped to the Eldorado urban extent.

Workflow and scripts
Script 1 – Environmental data preprocessing

Loads remote sensing rasters

Crops and masks layers using the Eldorado city polygon

Resamples all predictors to a common 10 m grid

Reprojects layers when needed

Outputs season-specific raster stacks used for modeling

Script 2 – Occurrence data preprocessing

Loads seasonal occurrence records

Removes records with missing or invalid coordinates

Removes duplicate geographic coordinates

Ensures a single presence point per raster cell

Extracts environmental values at presence locations

Outputs cleaned presence tables for each species and season

Script 3 – Variable selection

Performs multicollinearity analysis using:

Variance Inflation Factor (VIF)

Pearson correlation

Retains non-collinear predictors (VIF < 5)

Produces final predictor tables for model calibration

Script 4 – Model calibration and selection

Calibrates MaxEnt candidate models using:

Multiple feature class combinations

Regularization multipliers

Applies 10-fold random cross-validation

Evaluates models using ENMeval metrics

Selects the best model per species and season

Produces:

Model evaluation tables

Variable importance

Response curves

Script 5 – Model projection and binarization

Projects selected MaxEnt models to the study area

Uses cloglog output

Applies a 10th percentile training presence (p10) threshold

Generates continuous and binary suitability maps for each season and year

Script 6 – Seasonal aggregation and spatial analysis

Reprojects all maps to UTM Zone 21S

Averages continuous suitability maps by season across years

Computes seasonal differences (summer − winter)

Averages binary maps to represent the proportion of years classified as suitable

Calculates:

Percentage of suitable area per species and season

Suitable area in hectares within the Eldorado boundary

Outputs final maps and summary tables

Notes on interpretation

Averaged binary maps represent the temporal consistency of suitability, not strict presence/absence

Distance-to predictors capture spatial gradients of urbanization and vegetation influence

Nearest-neighbor interpolation is used for binary rasters; bilinear interpolation for continuous variables

Models are trained using presence and background data only (no true absences)

Output structure

raster/ – preprocessed environmental layers

data/ – cleaned occurrence tables and extracted predictors

best models/ – final MaxEnt predictions and binary maps

results/ – seasonal averages, difference maps, and area summaries

Reproducibility

All scripts are numbered and should be run sequentially.
File paths assume the repository root as the working directory.
