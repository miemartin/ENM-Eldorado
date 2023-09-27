require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, gtools, 
               corrplot, colourpicker, colorspace, usdm)

### Variables Selection

## Load data

pnts_albo <- read.csv("data/s2016_albo_clean.csv")
head(pnts_albo)

pnts_aeg <- read.csv("data/s2016_aeg_clean.csv")
head(pnts_aeg)

## VIF Analysis 

pnts_albo <- as.data.frame(pnts_albo)
vif.res <- vif(x = pnts_albo[,3:ncol(pnts_albo)])
vif.step <- vifstep(x = pnts_albo[,3:ncol(pnts_albo)], th = 5)
vrs <- vif.step@results$Variables %>% as.character()

pnts_aeg <- as.data.frame(pnts_aeg)
vif.res1 <- vif(x = pnts_aeg[,3:ncol(pnts_aeg)])
vif.step1 <- vifstep(x = pnts_aeg[,3:ncol(pnts_aeg)], th = 5) 
vrs1 <- vif.step1@results$Variables %>% as.character()


# Write the final table (vars selected)

pnts_albo <- dplyr::select(pnts_albo, x, y, all_of(vrs))
pnts_albo <- as_tibble(pnts_albo)
write.csv(pnts_albo, 'data/s2016_albo_clean_vars.csv', row.names = FALSE)

pnts_aeg <- dplyr::select(pnts_aeg, x, y, all_of(vrs1))
pnts_aeg <- as_tibble(pnts_aeg)
write.csv(pnts_aeg, 'data/s2016_aeg_clean_vars.csv', row.names = FALSE)

## Correlation

corr1 <- cor(pnts_albo[,3:ncol(pnts_albo)])
corr1
corrplot(corr1, method = 'square')

corr2 <- cor(pnts_aeg[,3:ncol(pnts_aeg)])
corr2
corrplot(corr2, method = 'square')
