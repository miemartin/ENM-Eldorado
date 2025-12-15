## =========================================================
## Script 1 – Environmental preprocessing
## Variable: Land Surface Temperature (LST)
## Season: Spring–Summer 2016
## =========================================================

require(pacman)
pacman::p_load(terra, sf)

## -----------------------------
## Parameters
## -----------------------------
season <- "s2016"
raster_dir <- "raster/"
shp_dir <- "shp/"

## -----------------------------
## Load rasters
## -----------------------------
water <- rast(paste0(raster_dir, season, "_water.tif"))
LST_30m <- rast(paste0(raster_dir, season, "_LST.tif"))

## -----------------------------
## Load and prepare study area
## -----------------------------
ED <- vect(paste0(shp_dir, "polygon_Eldorado.shp"))

## Ensure CRS consistency
ED <- project(ED, crs(LST_30m))

## -----------------------------
## Crop and mask LST
## -----------------------------
LST_crop <- crop(LST_30m, ED)
LST_mask <- mask(LST_crop, ED)

## -----------------------------
## Resample LST to 10 m
## -----------------------------
## Nearest neighbor used to preserve original LST values
LST_10m <- resample(LST_mask, water, method = "near")

## -----------------------------
## Save output
## -----------------------------
writeRaster(
  LST_10m,
  filename = paste0(raster_dir, season, "_LST_10m.tif"),
  overwrite = TRUE
)
