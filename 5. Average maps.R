library(terra)
library(raster)
library(rstudioapi)
library(dismo)

### AVERAGE PREDICTED MAPS

s2016_aeg <- raster(paste0("best models season/map_s2016_aeg.tif"))
s2017_aeg <- raster(paste0("best models season/map_s2017_aeg.tif"))
w2016_aeg <- raster(paste0("best models season/map_w2016_aeg.tif"))
w2017_aeg <- raster(paste0("best models season/map_w2017_aeg.tif"))

plot(s2016_aeg)
plot(s2017_aeg)

s2016_albo <- raster(paste0("best models season/map_s2016_albo.tif"))
s2017_albo <- raster(paste0("best models season/map_s2017_albo.tif"))
w2016_albo <- raster(paste0("best models season/map_w2016_albo.tif"))
w2017_albo <- raster(paste0("best models season/map_w2017_albo.tif"))

summer_aeg <- stack(s2016_aeg,s2017_aeg)
winter_aeg <- stack(w2016_aeg,w2017_aeg)

summer_albo <- stack(s2016_albo,s2017_albo)
winter_albo <- stack(w2016_albo,w2017_albo)

s_aeg <- mean(summer_aeg)
s_albo <- mean(summer_albo)
w_aeg <- mean(winter_aeg)
w_albo <- mean(winter_albo)

plot(s_aeg)

terra::writeRaster(x = s_aeg, filename = 'results/summer_aegypti.tif', overwrite=TRUE)
terra::writeRaster(x = s_albo, filename = 'results/summer_albopictus.tif', overwrite=TRUE)
terra::writeRaster(x = w_aeg, filename = 'results/winter_aegypti.tif', overwrite=TRUE)
terra::writeRaster(x = w_albo, filename = 'results/winter_albopictus.tif', overwrite=TRUE)


### AVERAGE BINARY MAPS 

s2016_aeg <- raster(paste0("best models season/map_s2016_aeg_p10.tif"))
s2017_aeg <- raster(paste0("best models season/map_s2017_aeg_p10.tif"))
w2016_aeg <- raster(paste0("best models season/map_w2016_aeg_p10.tif"))
w2017_aeg <- raster(paste0("best models season/map_w2017_aeg_p10.tif"))

plot(s2016_aeg)
plot(s2017_aeg)

s2016_albo <- raster(paste0("best models season/map_s2016_albo_p10.tif"))
s2017_albo <- raster(paste0("best models season/map_s2017_albo_p10.tif"))
w2016_albo <- raster(paste0("best models season/map_w2016_albo_p10.tif"))
w2017_albo <- raster(paste0("best models season/map_w2017_albo_p10.tif"))

summer_aeg <- stack(s2016_aeg,s2017_aeg)
winter_aeg <- stack(w2016_aeg,w2017_aeg)

summer_albo <- stack(s2016_albo,s2017_albo)
winter_albo <- stack(w2016_albo,w2017_albo)

s_aeg <- mean(summer_aeg)
s_albo <- mean(summer_albo)
w_aeg <- mean(winter_aeg)
w_albo <- mean(winter_albo)

plot(s_aeg)

terra::writeRaster(x = s_aeg, filename = 'results/s_b_aegypti.tif', overwrite=TRUE)
terra::writeRaster(x = s_albo, filename = 'results/s_b_albopictus.tif', overwrite=TRUE)
terra::writeRaster(x = w_aeg, filename = 'results/w_b_aegypti.tif', overwrite=TRUE)
terra::writeRaster(x = w_albo, filename = 'results/w_b_albopictus.tif', overwrite=TRUE)


### SUITABLE AREAS

# create empty matrix to be filled with suitable area values 
suitability <- matrix(data = NA, nrow = 4, ncol = 3)

# create object that lists the frequency of 0 and 1 pixels 
f <- freq(s_aeg)
# calculate percentage of area suitable by dividing number of 1 pixels by the number of 0 + 1 pixels
area1 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
# write results and label to the matrix 
suitability[1,1] <- paste0("Ae. aegypti")
suitability[1,2] <- paste0("summer")
suitability[1,3] <- area1

f <- freq(w_aeg)
area2 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[2,1] <- paste0("Ae. aegypti")
suitability[2,2] <- paste0("winter")
suitability[2,3] <- area2

f <- freq(s_albo)
area3 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[3,1] <- paste0("Ae. albopictus")
suitability[3,2] <- paste0("summer")
suitability[3,3] <- area3

f <- freq(w_albo)
area4 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[4,1] <- paste0("Ae. albopictus")
suitability[4,2] <- paste0("winter")
suitability[4,3] <- area4

# format table and write .csv to the BinaryMaps folder 
colnames(suitability) <- c("Specie","Season","SuitableArea")
write.csv(suitability, file = "SuitableArea.csv")
