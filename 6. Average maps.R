## =========================================================
## Script 6 – Seasonal aggregation and spatial analysis
## =========================================================

library(terra)
library(sf)
library(tidyverse)

## -----------------------------
## CRS definition
## -----------------------------
crs_utm21s <- "EPSG:32721"  # UTM zone 21S, WGS84

## -----------------------------
## Load continuous predictions
## -----------------------------
s2016_aeg <- rast("best models/map_s2016_aegypti.tif")
s2017_aeg <- rast("best models/map_s2017_aegypti.tif")
w2016_aeg <- rast("best models/map_w2016_aegypti.tif")
w2017_aeg <- rast("best models/map_w2017_aegypti.tif")

s2016_albo <- rast("best models/map_s2016_albo.tif")
s2017_albo <- rast("best models/map_s2017_albo.tif")
w2016_albo <- rast("best models/map_w2016_albo.tif")
w2017_albo <- rast("best models/map_w2017_albo.tif")

## Reproject (continuous → bilinear)
s2016_aeg <- project(s2016_aeg, crs_utm21s, method = "bilinear")
s2017_aeg <- project(s2017_aeg, crs_utm21s, method = "bilinear")
w2016_aeg <- project(w2016_aeg, crs_utm21s, method = "bilinear")
w2017_aeg <- project(w2017_aeg, crs_utm21s, method = "bilinear")

s2016_albo <- project(s2016_albo, crs_utm21s, method = "bilinear")
s2017_albo <- project(s2017_albo, crs_utm21s, method = "bilinear")
w2016_albo <- project(w2016_albo, crs_utm21s, method = "bilinear")
w2017_albo <- project(w2017_albo, crs_utm21s, method = "bilinear")

## -----------------------------
## Seasonal averages (continuous)
## -----------------------------
summer_aeg <- mean(c(s2016_aeg, s2017_aeg))
winter_aeg <- mean(c(w2016_aeg, w2017_aeg))

summer_albo <- mean(c(s2016_albo, s2017_albo))
winter_albo <- mean(c(w2016_albo, w2017_albo))

writeRaster(summer_aeg, "results/summer_aegypti.tif", overwrite = TRUE)
writeRaster(winter_aeg, "results/winter_aegypti.tif", overwrite = TRUE)
writeRaster(summer_albo, "results/summer_albopictus.tif", overwrite = TRUE)
writeRaster(winter_albo, "results/winter_albopictus.tif", overwrite = TRUE)

## -----------------------------
## Seasonal differences
## -----------------------------
writeRaster(summer_aeg - winter_aeg, "results/dif_sw_aegypti.tif", overwrite = TRUE)
writeRaster(winter_aeg - summer_aeg, "results/dif_ws_aegypti.tif", overwrite = TRUE)

writeRaster(summer_albo - winter_albo, "results/dif_sw_albopictus.tif", overwrite = TRUE)
writeRaster(winter_albo - summer_albo, "results/dif_ws_albopictus.tif", overwrite = TRUE)

## -----------------------------
## Load binary maps (p10)
## -----------------------------
s2016_aeg_b <- rast("best models/map_s2016_aeg_p10.tif")
s2017_aeg_b <- rast("best models/map_s2017_aeg_p10.tif")
w2016_aeg_b <- rast("best models/map_w2016_aeg_p10.tif")
w2017_aeg_b <- rast("best models/map_w2017_aeg_p10.tif")

s2016_albo_b <- rast("best models/map_s2016_albo_p10.tif")
s2017_albo_b <- rast("best models/map_s2017_albo_p10.tif")
w2016_albo_b <- rast("best models/map_w2016_albo_p10.tif")
w2017_albo_b <- rast("best models/map_w2017_albo_p10.tif")

## Reproject (binary → nearest neighbor)
s2016_aeg_b <- project(s2016_aeg_b, crs_utm21s, method = "near")
s2017_aeg_b <- project(s2017_aeg_b, crs_utm21s, method = "near")
w2016_aeg_b <- project(w2016_aeg_b, crs_utm21s, method = "near")
w2017_aeg_b <- project(w2017_aeg_b, crs_utm21s, method = "near")

s2016_albo_b <- project(s2016_albo_b, crs_utm21s, method = "near")
s2017_albo_b <- project(s2017_albo_b, crs_utm21s, method = "near")
w2016_albo_b <- project(w2016_albo_b, crs_utm21s, method = "near")
w2017_albo_b <- project(w2017_albo_b, crs_utm21s, method = "near")

## -----------------------------
## Average binary maps
## (proportion of years suitable)
## -----------------------------
s_aeg_b <- mean(c(s2016_aeg_b, s2017_aeg_b))
w_aeg_b <- mean(c(w2016_aeg_b, w2017_aeg_b))

s_albo_b <- mean(c(s2016_albo_b, s2017_albo_b))
w_albo_b <- mean(c(w2016_albo_b, w2017_albo_b))

writeRaster(s_aeg_b, "results/s_b_aegypti.tif", overwrite = TRUE)
writeRaster(w_aeg_b, "results/w_b_aegypti.tif", overwrite = TRUE)
writeRaster(s_albo_b, "results/s_b_albopictus.tif", overwrite = TRUE)
writeRaster(w_albo_b, "results/w_b_albopictus.tif", overwrite = TRUE)

## -----------------------------
## Suitable area percentages
## -----------------------------
calc_area <- function(r) {
  f <- freq(r, useNA = "no")
  round(100 * f$value[f$value == 1] / sum(f$count), 2)
}

suitability <- tibble(
  Species = c("Ae. aegypti","Ae. aegypti","Ae. albopictus","Ae. albopictus"),
  Season  = c("summer","winter","summer","winter"),
  SuitableArea = c(
    calc_area(s_aeg_b),
    calc_area(w_aeg_b),
    calc_area(s_albo_b),
    calc_area(w_albo_b)
  )
)

write.csv(suitability, "SuitableArea.csv", row.names = FALSE)

## -----------------------------
## Area in hectares (example)
## -----------------------------
ED <- vect(st_transform(read_sf("shp/polygon_Eldorado.shp"), crs_utm21s))

pol <- as.polygons(s_aeg_b, dissolve = TRUE)
intr <- intersect(pol, ED)
intr$ha <- round(expanse(intr, unit = "ha"), 1)
