## =========================================================
## Script 3 – Multicollinearity assessment (VIF analysis)
## =========================================================

require(pacman)
pacman::p_load(tidyverse, usdm, corrplot)

## -----------------------------
## Parameters
## -----------------------------
season <- "s2016"
data_dir <- "data/"
vif_threshold <- 5

## -----------------------------
## Load data
## -----------------------------
albo <- read.csv(paste0(data_dir, season, "_albo_clean.csv"))
aeg  <- read.csv(paste0(data_dir, season, "_aeg_clean.csv"))

## Predictor columns (exclude coordinates)
predictor_cols <- setdiff(names(albo), c("x", "y"))

## -----------------------------
## Function for VIF selection
## -----------------------------
run_vif <- function(df, predictors, threshold) {

  df_sel <- df %>%
    dplyr::select(all_of(predictors)) %>%
    drop_na()

  vif_step <- vifstep(df_sel, th = threshold)

  selected_vars <- vif_step@results$Variables %>% as.character()

  return(selected_vars)
}

## -----------------------------
## Run VIF analysis
## -----------------------------
vars_albo <- run_vif(albo, predictor_cols, vif_threshold)
vars_aeg  <- run_vif(aeg,  predictor_cols, vif_threshold)

## -----------------------------
## Save reduced datasets
## -----------------------------
albo_final <- albo %>%
  dplyr::select(x, y, all_of(vars_albo))

aeg_final <- aeg %>%
  dplyr::select(x, y, all_of(vars_aeg))

write.csv(albo_final,
          paste0(data_dir, season, "_albo_clean_vars.csv"),
          row.names = FALSE)

write.csv(aeg_final,
          paste0(data_dir, season, "_aeg_clean_vars.csv"),
          row.names = FALSE)

## -----------------------------
## Save selected variables
## -----------------------------
write.csv(data.frame(variable = vars_albo),
          paste0(data_dir, season, "_albo_selected_variables.csv"),
          row.names = FALSE)

write.csv(data.frame(variable = vars_aeg),
          paste0(data_dir, season, "_aeg_selected_variables.csv"),
          row.names = FALSE)

## -----------------------------
## Correlation plots (diagnostic)
## -----------------------------
corr_albo <- cor(albo_final[, -c(1,2)])
corr_aeg  <- cor(aeg_final[, -c(1,2)])

corrplot(corr_albo, method = "square")
corrplot(corr_aeg,  method = "square")
