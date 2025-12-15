## =========================================================
## Script 2 – Occurrence data processing and variable extraction
## Season: Spring–Summer 2016
## =========================================================

require(pacman)
pacman::p_load(terra, tidyverse)

## -----------------------------
## Parameters
## -----------------------------
season <- "s2016"
raster_dir <- "raster/"
data_dir   <- "data/"

## -----------------------------
## Load environmental variables
## -----------------------------
vars <- rast(c(
  paste0(raster_dir, season, "_water.tif"),
  paste0(raster_dir, season, "_urban.tif"),
  paste0(raster_dir, season, "_bare.tif"),
  paste0(raster_dir, season, "_tree.tif"),
  paste0(raster_dir, season, "_herb.tif"),
  paste0(raster_dir, season, "_dwater.tif"),
  paste0(raster_dir, season, "_dimp.tif"),
  paste0(raster_dir, season, "_dveg.tif"),
  paste0(raster_dir, season, "_NDVI.tif"),
  paste0(raster_dir, season, "_NDBI.tif"),
  paste0(raster_dir, season, "_NDWI.tif"),
  paste0(raster_dir, season, "_LST.tif")
))

names(vars) <- c("water","urban","soil","tree","herb",
                 "d_water","d_imp","d_veg",
                 "NDVI","NDBI","NDWI","LST")

## -----------------------------
## Load occurrence data
## -----------------------------
occ_albo_raw <- read.csv(paste0(data_dir, season, "_albo.csv"), sep = ";")
occ_aeg_raw  <- read.csv(paste0(data_dir, season, "_aeg.csv"),  sep = ";")

## -----------------------------
## Clean occurrence data
## -----------------------------
clean_occ <- function(df) {
  df %>%
    filter(!is.na(latitude), !is.na(longitude)) %>%
    distinct(latitude, longitude, .keep_all = TRUE)
}

occ_albo <- clean_occ(occ_albo_raw)
occ_aeg  <- clean_occ(occ_aeg_raw)

## -----------------------------
## Convert to spatial points
## -----------------------------
pts_albo <- vect(occ_albo, geom = c("longitude", "latitude"), crs = crs(vars))
pts_aeg  <- vect(occ_aeg,  geom = c("longitude", "latitude"), crs = crs(vars))

## -----------------------------
## One occurrence per raster cell
## -----------------------------
cells_albo <- cellFromXY(vars[[1]], geom(pts_albo)[,c("x","y")])
cells_aeg  <- cellFromXY(vars[[1]], geom(pts_aeg)[,c("x","y")])

pts_albo <- pts_albo[!duplicated(cells_albo), ]
pts_aeg  <- pts_aeg[!duplicated(cells_aeg),  ]

## -----------------------------
## Extract environmental values
## -----------------------------
albo_vals <- terra::extract(vars, pts_albo, xy = TRUE)
aeg_vals  <- terra::extract(vars, pts_aeg,  xy = TRUE)

albo_vals <- as_tibble(albo_vals)
aeg_vals  <- as_tibble(aeg_vals)

## -----------------------------
## Save outputs
## -----------------------------
write.csv(albo_vals, paste0(data_dir, season, "_albo_clean.csv"), row.names = FALSE)
write.csv(aeg_vals,  paste0(data_dir, season, "_aeg_clean.csv"),  row.names = FALSE)
