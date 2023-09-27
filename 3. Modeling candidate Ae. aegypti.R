library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(raster)
library(rgdal)
library(ENMeval)
library(wallace)
library(sf)
library(stats)

### Analysis for Aedes aegypti (Aa)

## Load occurrence data

occs_path <- "C:/Users/Giga/Desktop/aedes potential maps/Eldorado/season models/season models/data/"
occs_path <- file.path(occs_path, "wallace_aeg_s2017.csv")

userOccs_Aa <- occs_userOccs(
  txtPath = occs_path, 
  txtName = "wallace_aeg_s2017.csv", 
  txtSep = ",", 
  txtDec = ".")
occs_Aa <- userOccs_Aa$Aedes_aegypti$cleaned


## Load remote sensing data

dir_envs_Aa <- "C:/Users/Giga/Desktop/aedes potential maps/Eldorado/season models/season models/raster"
envs_path <- file.path(dir_envs_Aa, c("s2017_water.tif","s2017_bare.tif","s2017_herb.tif","s2017_dwater.tif","s2017_dimp.tif", "s2017_dveg.tif", "s2017_NDVI.tif","s2017_NDWI.tif", "s2017_LST.tif"))

envs_Aa <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c("s2017_water.tif","s2017_bare.tif","s2017_herb.tif","s2017_dwater.tif","s2017_dimp.tif", "s2017_dveg.tif", "s2017_NDVI.tif","s2017_NDWI.tif", "s2017_LST.tif"),
  doBrick = FALSE)

occs_xy_Aa <- occs_Aa[c('longitude', 'latitude')]
occs_vals_Aa <- as.data.frame(raster::extract(envs_Aa, occs_xy_Aa, cellnumbers = TRUE))

# add columns for env variable values for each occurrence record
occs_Aa <- cbind(occs_Aa, occs_vals_Aa)


### Process remote sensing data
#Sampling of 100 background points and corresponding environmental data using a “point buffers” method with a 2 degree buffer.

# Generate background extent 
bgExt_Aa <- penvs_bgExtent(
  occs = occs_Aa,
  bgSel = "point buffers",
  bgBuf = 0.008)

# Mask environmental data to provided extent
bgMask_Aa <- penvs_bgMask(
  occs = occs_Aa,
  envs = envs_Aa,
  bgExt = bgExt_Aa)

# Sample background points from the provided area
bgSample_Aa <- penvs_bgSample(
  occs = occs_Aa,
  bgMask =  bgMask_Aa,
  bgPtsNum = 10000)

# Save dataframe to shp
WGScoor<-  bgExt_Aa
coordinates(WGScoor)=~longitude+latitude
proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
bgcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
raster::shapefile(bgcoor, "shp/s2017_Marea_aegypti.shp")

WGScoor<-  bgSample_Aa
coordinates(WGScoor)=~longitude+latitude
proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
bgcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
raster::shapefile(bgcoor, "shp/s2017_bg_aegypti.shp")

WGScoor<-  occs_Aa
coordinates(WGScoor)=~longitude+latitude
proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
coor<-spTransform(WGScoor,CRS("+proj=longlat"))
raster::shapefile(coor, "shp/s2017_occ_aeg.shp")


# Extract values of environmental layers for each background point
bgEnvsVals_Aa <- as.data.frame(raster::extract(bgMask_Aa,  bgSample_Aa))

# Add extracted values to background points table
bgEnvsVals_Aa <- cbind(scientific_name = paste0("bg_", "Aedes aegypti"), bgSample_Aa,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Aa)


### Partition occurrence data
#Partition occurrences and background points for model training and validation using random k-fold, a non-spatial partition method.

groups_Aa <- part_partitionOccs(
  occs = occs_Aa ,
  bg =  bgSample_Aa, 
  method = "rand",
  kfolds = 2) 


### Build and Evaluate Niche Model
# Generating a species distribution model using the maxent.jar algorithm as implemented in ENMeval V2.0 (with clamping = TRUE). 
# For tuning using L, LQ, LQH, LQHP feature classes.

## Run maxent model for the selected species

# regularization multipliers in the 0.5-1 range increasing by 0.1.
model_Aa <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(0.1, 1), 
  rmsStep =  0.1,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  parallel = FALSE)

evalTbl <- model_Aa@results 
write.csv(evalTbl, "s2017_candidate_models_aegypti1.csv", row.names=FALSE)

# regularization multipliers in the 2-6 range increasing by 1. 

model_Aa2 <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(2, 6), 
  rmsStep =  1,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  parallel = FALSE)

evalTbl2 <- model_Aa2@results 
write.csv(evalTbl2, "s2017_candidate_models_aegypti2.csv", row.names=FALSE)

# regularization multipliers 8 and 10

model_Aa3 <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(8, 10), 
  rmsStep =  2,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  parallel = FALSE)

evalTbl3 <- model_Aa3@results 
write.csv(evalTbl3, "s2017_candidate_models_aegypti3.csv", row.names=FALSE)
