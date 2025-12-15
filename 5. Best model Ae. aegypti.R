## =========================================================
## Script 5 – Final MaxEnt model (Ae. aegypti)
## =========================================================

require(pacman)
pacman::p_load(ENMeval, wallace, dismo, terra, sf)

## -----------------------------
## Notes on inputs
## -----------------------------
## The following objects are assumed to be preprocessed and
## loaded from previous scripts:
## - occs_Aa        : occurrence records (data.frame)
## - occs_xy_Aa     : coordinates of occurrences
## - bgSample_Aa    : background points
## - bgEnvsVals_Aa  : background environmental values
## - bgMask_Aa      : raster mask for prediction
## - envs_Aa        : environmental rasters for transfer

set.seed(123)

## -----------------------------
## Partition occurrence data
## -----------------------------
groups_Aa <- part_partitionOccs(
  occs = occs_Aa,
  bg   = bgSample_Aa,
  method = "rand",
  kfolds = 10
)

## -----------------------------
## Run best-selected model
## -----------------------------
model_Aa_best <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa,
  bgMsk = bgMask_Aa,
  rms = 5,
  fcs = "LQH",
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  catEnvs = c("s2017_water","s2017_bare","s2017_herb"),
  parallel = FALSE
)

## -----------------------------
## Save evaluation outputs
## -----------------------------
write.csv(model_Aa_best@results,
          "best_s2016_aegypti.csv",
          row.names = FALSE)

write.csv(model_Aa_best@variable.importance,
          "imp_s2017_aeg.csv",
          row.names = FALSE)

## Response curves
response(model_Aa_best@models[["fc.LQH_rm.5"]])

## -----------------------------
## Transfer model to study area
## -----------------------------
xfer_extent <- vect("shp/polygon_Eldorado.shp")

xfer_out <- xfer_area(
  evalOut = model_Aa_best,
  curModel = "fc.LQH_rm.5",
  envs = envs_Aa,
  outputType = "cloglog",
  alg = "maxent.jar",
  clamp = TRUE,
  xfExt = xfer_extent
)

writeRaster(xfer_out$xferArea,
            "map_s2017_aegypti.tif",
            overwrite = TRUE)

## -----------------------------
## Binary map (p10 threshold)
## -----------------------------
m_Aa <- model_Aa_best@models[["fc.LQH_rm.5"]]

pred <- predict(
  m_Aa,
  bgMask_Aa,
  args = c("outputformat=cloglog", "doclamp=true"),
  na.rm = TRUE
)

occ_vals <- raster::extract(pred, occs_xy_Aa)

thres_Aa <- quantile(occ_vals, probs = 0.1, na.rm = TRUE)

pred_bin <- pred > thres_Aa

writeRaster(pred_bin,
            "map_s2017_aegypti_p10.tif",
            overwrite = TRUE)

## -----------------------------
## Binary transfer
## -----------------------------
xfer_bin <- xfer_out$xferArea > thres_Aa

writeRaster(xfer_bin,
            "map_s2017_aegypti_xfer_p10.tif",
            overwrite = TRUE)
