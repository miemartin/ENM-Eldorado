## =========================================================
## Script 4 – Seasonal averaging, differences, and area estimates
## =========================================================

require(pacman)
pacman::p_load(terra, sf, tidyverse)

## -----------------------------
## Parameters
## -----------------------------
crs_utm <- "EPSG:32721"  # UTM zone 21S
model_dir <- "best models/"
res_dir <- "results/"
shp_dir <- "shp/"

## -----------------------------
## Helper function
## -----------------------------
load_and_project <- function(files, crs_out) {
  rast(files) |> project(crs_out)
}

## -----------------------------
## Load continuous suitability maps
## -----------------------------
aeg_summer <- load_and_project(
  c("map_s2016_aegypti.tif", "map_s2017_aegypti.tif") |> 
    file.path(model_dir, .),
  crs_utm
)

aeg_winter <- load_and_project(
  c("map_w2016_aegypti.tif", "map_w2017_aegypti.tif") |> 
    file.path(model_dir, .),
  crs_utm
)

albo_summer <- load_and_project(
  c("map_s2016_albo.tif", "map_s2017_albo.tif") |> 
    file.path(model_dir, .),
  crs_utm
)

albo_winter <- load_and_project(
  c("map_w2016_albo.tif", "map_w2017_albo.tif") |> 
    file.path(model_dir, .),
  crs_utm
)

## -----------------------------
## Seasonal mean suitability
## -----------------------------
s_aeg  <- mean(aeg_summer)
w_aeg  <- mean(aeg_winter)
s_albo <- mean(albo_summer)
w_albo <- mean(albo_winter)

writeRaster(s_aeg,  file.path(res_dir, "summer_aegypti.tif"), overwrite = TRUE)
writeRaster(w_aeg,  file.path(res_dir, "winter_aegypti.tif"), overwrite = TRUE)
writeRaster(s_albo, file.path(res_dir, "summer_albopictus.tif"), overwrite = TRUE)
writeRaster(w_albo, file.path(res_dir, "winter_albopictus.tif"), overwrite = TRUE)

## -----------------------------
## Seasonal differences
## -----------------------------
writeRaster(s_aeg - w_aeg,  file.path(res_dir, "dif_sw_aegypti.tif"), overwrite = TRUE)
writeRaster(w_aeg - s_aeg,  file.path(res_dir, "dif_ws_aegypti.tif"), overwrite = TRUE)
writeRaster(s_albo - w_albo,file.path(res_dir, "dif_sw_albopictus.tif"), overwrite = TRUE)
writeRaster(w_albo - s_albo,file.path(res_dir, "dif_ws_albopictus.tif"), overwrite = TRUE)

## -----------------------------
## Load binary (p10) maps
## -----------------------------
load_binary <- function(files) {
  rast(files |> file.path(model_dir, .))
}

aeg_s_bin <- mean(load_binary(c("map_s2016_aeg_p10.tif","map_s2017_aeg_p10.tif")))
aeg_w_bin <- mean(load_binary(c("map_w2016_aeg_p10.tif","map_w2017_aeg_p10.tif")))
albo_s_bin <- mean(load_binary(c("map_s2016_albo_p10.tif","map_s2017_albo_p10.tif")))
albo_w_bin <- mean(load_binary(c("map_w2016_albo_p10.tif","map_w2017_albo_p10.tif")))

writeRaster(aeg_s_bin,  file.path(res_dir,"s_b_aegypti.tif"), overwrite=TRUE)
writeRaster(aeg_w_bin,  file.path(res_dir,"w_b_aegypti.tif"), overwrite=TRUE)
writeRaster(albo_s_bin, file.path(res_dir,"s_b_albopictus.tif"), overwrite=TRUE)
writeRaster(albo_w_bin, file.path(res_dir,"w_b_albopictus.tif"), overwrite=TRUE)

## -----------------------------
## Suitable area percentages
## -----------------------------
calc_area <- function(r) {
  f <- freq(r, useNA = "no")
  round(100 * f$value[f$value==1] / sum(f$count), 2)
}

suitability <- tibble(
  Species = c("Ae. aegypti","Ae. aegypti","Ae. albopictus","Ae. albopictus"),
  Season  = c("summer","winter","summer","winter"),
  SuitableArea = c(
    calc_area(aeg_s_bin),
    calc_area(aeg_w_bin),
    calc_area(albo_s_bin),
    calc_area(albo_w_bin)
  )
)

write.csv(suitability, "SuitableArea.csv", row.names = FALSE)

## -----------------------------
## Area in hectares (optional)
## -----------------------------
ED <- vect(file.path(shp_dir,"polygon_Eldorado.shp")) |> project(crs_utm)

ha_map <- as.polygons(aeg_s_bin) |> intersect(ED)
ha_map$ha <- round(expanse(ha_map) / 10000, 1)
