require(pacman)
pacman::p_load(geodata, terra, glue, rgeos, gtools, ggspatial, sf, 
               colourpicker, tidyverse, fs, rnaturalearthdata, 
               rnaturalearth, glue, outliers, stringr, raster)

### Load remote sensing data (spring-summer 2016)

water <- raster(paste0("raster/s2016_water.tif"))
urban <- raster(paste0("raster/s2016_urban.tif"))
soil <- raster(paste0("raster/s2016_bare.tif"))
tree <- raster(paste0("raster/s2016_tree.tif"))
herb <- raster(paste0("raster/s2016_herb.tif"))
d_water <- raster(paste0("raster/s2016_dwater.tif"))
d_imp <- raster(paste0("raster/s2016_dimp.tif"))
d_veg <- raster(paste0("raster/s2016_dveg.tif"))
NDVI <- raster(paste0("raster/s2016_NDVI.tif"))
NDBI <- raster(paste0("raster/s2016_NDBI.tif"))
NDWI <- raster(paste0("raster/s2016_NDWI.tif"))
LST <- raster(paste0("raster/s2016_LST.tif"))

# stacking the variables 
var <- raster::stack(water,urban,soil,tree,herb,d_water,d_imp,d_veg,NDVI,NDBI,NDWI,LST)


### Load occurrence data (spring-summer 2016)

occ_albo_raw <- read.csv("data/s2016_albo.csv", sep=";")
occ_aeg_raw <- read.csv("data/s2016_aeg.csv", sep=";")


### Clean occurrence data

# remove erroneous coordinates, where either the latitude or longitude is missing
occ_albo_clean <- subset(occ_albo_raw,(!is.na(latitude))&(!is.na(longitude)))
cat(nrow(occ_albo_raw)-nrow(occ_albo_clean), "records are removed") 

occ_aeg_clean <- subset(occ_aeg_raw,(!is.na(latitude))&(!is.na(longitude)))
cat(nrow(occ_aeg_raw)-nrow(occ_aeg_clean), "records are removed") 

# remove duplicated data based on latitude and longitude
dups_albo <- duplicated(occ_albo_clean[c("latitude","longitude")])
occ_albo <- occ_albo_clean[!dups_albo,]
cat(nrow(occ_albo_clean)-nrow(occ_albo), "records are removed") 

dups_aeg <- duplicated(occ_aeg_clean[c("latitude","longitude")])
occ_aeg <- occ_aeg_clean[!dups_aeg,]
cat(nrow(occ_aeg_clean)-nrow(occ_aeg), "records are removed") 

# Only one occurrence point per pixel

coordinates(occ_albo) <- ~ longitude + latitude
cells1 <- cellFromXY(var[[1]],occ_albo)
dups1 <- duplicated(cells1)
occ_final1 <- occ_albo[!dups1,]
cat(nrow(occ_albo)-nrow(occ_final1), "records are removed") 

coordinates(occ_aeg) <- ~ longitude + latitude
cells2 <- cellFromXY(var[[1]],occ_aeg)
dups2 <- duplicated(cells2)
occ_final2 <- occ_aeg[!dups2,]
cat(nrow(occ_aeg)-nrow(occ_final2), "records are removed")

# Write the result as table
albo <- write.csv(occ_final1, 'data/albopictus_final.csv', row.names = FALSE) 
aeg <- write.csv(occ_final2, 'data/aegypti_final.csv', row.names = FALSE) 

# Extract the values for the presences 
vles1 <- terra::extract(var, albo[,2:3])
vles1 <- as_tibble(cbind(albo, vles1))

vles2 <- terra::extract(var, aeg[,2:3])
vles2 <- as_tibble(cbind(aeg, vles2))

# Write the result (presences and remote sensing data)
write.csv(vles1, 'data/s2016_albo_vls.csv', row.names = FALSE)
write.csv(vles2, 'data/s2016_aeg_vls.csv', row.names = FALSE)

# Add remote sensing data values
pnts1 <- cbind(vles1[,2:3], raster::extract(var, vles1[,2:3]))
names(pnts1) <- c("x", "y", "water","urban","soil","tree","herb","d_water","d_imp","d_veg","NDVI","NDBI","NDWI","LST")
pnts1 <- as_tibble(pnts1)
write.csv(pnts1, 'data/s2016_albo_clean.csv', row.names = FALSE)


pnts2 <- cbind(vles2[,2:3], raster::extract(var, vles2[,2:3]))
names(pnts2) <- c("x", "y", "water","urban","soil","tree","herb","d_water","d_imp","d_veg","NDVI","NDBI","NDWI","LST")
pnts2 <- as_tibble(pnts2)
write.csv(pnts2, 'data/s2016_aeg_clean.csv', row.names = FALSE)