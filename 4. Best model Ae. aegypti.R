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

### Partition occurrence data

## Selected model: fc.LQHP_rm.4

groups_Aa2 <- part_partitionOccs(
  occs = occs_Aa ,
  bg =  bgSample_Aa, 
  method = "rand",
  kfolds = 10) 

## Run best model

model_Aa_best <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa2, 
  bgMsk = bgMask_Aa,
  rms = c(4, 5), 
  rmsStep =  1.5,
  fcs = c("LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  parallel = FALSE)

evalTbl_b <- model_Aa_best@results
write.csv(evalTbl_b, "best_s2016_aegypti.csv", row.names=FALSE)
evalMods <- model_Aa_best@models
evalPred <- model_Aa_best@predictions
evalVar <- model_Aa_best@variable.importance
write.csv(evalVar, "imp_s2016_aeg.csv", row.names=FALSE)

# view response curves for environmental variables
response(evalMods[["fc.LQHP_rm.4"]])


### Transfer model

#Transfering the model to a new user provided area

xfer_userExt_Aa <- read_sf("shp/polygon_Eldorado.shp")

# Generate a transfer of the model to the desired area
xfer_area_Aa <- xfer_area(
  evalOut = model_Aa_best,
  curModel = "fc.LQHP_rm.4",
  envs = envs_Aa , 
  outputType = "cloglog",
  alg = "maxent.jar",
  clamp = TRUE,
  xfExt = xfer_userExt_Aa) 

#store the cropped transfer variables
xferExt_Aa <- xfer_area_Aa$xferExt

plot(xfer_area_Aa$xferArea)
writeRaster(xfer_area_Aa$xferArea, "map_s2016_aegypti.tif")


### Binary map

#Generate a map of the MaxEnt generated model with a “p10” threshold rule.
#10 percentile training presence: 

# Select current model and obtain raster prediction
m_Aa <- model_Aa_best@models[["fc.LQHP_rm.4"]] 
predSel_Aa <- dismo::predict(m_Aa, bgMask_Aa,
                             args = c(paste0("outputformat=", "cloglog"), 
                                      paste0("doclamp=", tolower(as.character(TRUE)))), 
                             na.rm = TRUE)

# Determine the threshold based on the current prediction
occPredVals_Aa <- raster::extract(predSel_Aa, occs_xy_Aa)

# Define probability of quantile based on selected threshold
thresProb_Aa <- switch("p10", "mtp" = 0, "p10" = 0.1, "qtp" = 0)

# Define threshold value
thres_Aa <- stats::quantile(occPredVals_Aa, probs = thresProb_Aa, na.rm = TRUE)

# Applied selected threshold
predSel_Aa <- predSel_Aa > thres_Aa


## Transfer model
#Transfering the model to a new user drawn area using a “p10” threshold

# Add threshold if specified 
xfer_area_Aa <- xfer_area_Aa$xferArea > thres_Aa

# Plot and save the transfer raster
xferExt_Aa <- xfer_area_Aa
plot(xfer_area_Aa)
writeRaster(xfer_area_Aa, "map__s2016aeg_p10.tif")
