# R scripts for the course "Species distribution and ecological niche modelling in R"
# A. Marcia Barbosa (https://modtools.wordpress.com)


# don't forget to DO THIS FIRST:
# Session -> Set Working Directory -> To Source File Location


# LOAD PACKAGES ####

source("./scripts/00_packages.R")

library(dismo)
library(plotmo)
library(maxnet)
library(gam)
library(randomForest)
library(gbm)
library(maps)
library(fuzzySim)
library(modEvA)
library(sdmpredictors)
library(ecospat)
library(tidyverse)


# IMPORT DATA AND DEFINE SOME SETTINGS ####

# define the target species:
myspecies <- c("Linum bienne")

# import models saved in previous practical:
load(paste0("./data/niche_model/4modelling/", str_remove(myspecies, "Linum "), ".RData"))

# define colours and colour breaks for the maps:
clrs <- hcl.colors(10)
brks <- seq(0, 1, by = 0.1)


# PLOT VARIABLE RESPONSE CURVES ####

# before you project your models, you may want to take a look at the response curves, to see if they make ecological sense
# but mind that individual response curves may be affected by interactions with other variables

?response  # for 'dismo' models
dismo::response(bioclim_mod)
response(domain_mod)
response(maxent_mod)

?plot.maxnet  # for 'maxnet' models
plot(maxnet_mod)

?plotmo  # for other models
plotmo(glm_mod, trace = 1, all1 = TRUE)  # see the "plots" pane
plotmo(glm_mod, trace = 1, all1 = TRUE, all2 = TRUE)  # with all interactions
plotmo(gam_mod, trace = 1, all1 = TRUE)
plotmo(gam_mod, trace = 1, all1 = TRUE, all2 = TRUE)
plotmo(rf_mod, trace = 1, type = "prob", all1 = TRUE)
plotmo(rf_mod, trace = 1, type = "prob", all1 = TRUE, all2 = TRUE)
plotmo(gbm_mod, trace = 1, all1 = TRUE)
plotmo(gbm_mod, trace = 1, all1 = TRUE, all2 = TRUE)


# PROJECT THE MODELS TO A DIFFERENT REGION ####

# import the occurrences data:
occurrences <- read.csv(paste0("./data/niche_model/species_occurrence/", str_remove(myspecies, "Linum "), "_gbif_aug2022_clean.csv"))

## Also clip occurrences to geographical area of interest
my_window <- c(-22, 60, 25, 65)

occurrences <- occurrences %>% 
  filter(decimalLongitude >= my_window[1] & decimalLongitude <= my_window[2]) %>% #longitude
  droplevels() %>%
  filter(decimalLatitude >= my_window[3] & decimalLatitude <= my_window[4]) %>% #latitude
  droplevels() %>%
  as.data.frame()

# convert to a spatial object:
occurrences_spatial <- occurrences
names(occurrences_spatial)
coordinates(occurrences_spatial) <- occurrences[ , c("decimalLongitude", "decimalLatitude")]
crs(occurrences_spatial) <- "+proj=longlat"
par(mfrow = c(1, 1))
plot(occurrences_spatial)

# import the countries map:
countries <- readOGR("./data/niche_model/countries/world_countries.shp")
plot(countries, xlim = range(occurrences$decimalLongitude), ylim = range(occurrences$decimalLatitude))
map.axes()
points(occurrences_spatial, col = "blue")

# import the complete (global coverage) layers:
layers <- stack(list.files("./data/niche_model/climate/", pattern = "\\.tif$", full.names = TRUE))
plot(layers)

# now choose the projection extent, i.e. the region where you want to project your model predictions
# you can use e.g. the wider area that contains species occurrence points (including those that were left out of the modelling region):
#proj_extent <- extent(occurrences_spatial)
# or you can choose a specific spatial window as I did using my_window // (#DIDN't do it, but might use a buffer or draw a polygon for my species from Diederichsen paper)
proj_extent <- extent(c(-22, 60, 25, 65))

# crop the global layers to the projection extent:
layers_proj <- crop(layers, proj_extent)
plot(layers_proj, col = clrs)

# predict with each model to the projection layers:
bioclim_proj <- predict(layers_proj, bioclim_mod)
domain_proj <- predict(layers_proj[[names(domain_mod@presence)]], domain_mod)  # takes time!
maxent_proj <- predict(layers_proj, maxent_mod)
maxnet_proj <- predict(layers_proj, maxnet_mod, type = "cloglog")
glm_proj <- predict(layers_proj, glm_mod, type = "response")
gam_proj <- predict(layers_proj, gam_mod, type = "response")
rf_proj <- 1 - predict(layers_proj, rf_mod, type = "prob") # "1-" because I need occurrence probs
gbm_proj <- predict(layers_proj, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_fav <- Fav(pred = glm_proj, sample.preval = prevalence(glm_mod$y))
gam_proj_fav <- Fav(pred = gam_proj, sample.preval = prevalence(gam_mod$y))
rf_proj_fav <- Fav(pred = rf_proj, sample.preval = prevalence(rf_mod$y))
gbm_proj_fav <- Fav(pred = gbm_proj, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):

par(mfrow = c(4, 2))

plot(bioclim_proj, col = clrs, breaks = brks, main = "Bioclim")
plot(domain_proj, col = clrs, breaks = brks, main = "Domain") #does not work?
plot(maxent_proj, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_fav, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_fav, col = clrs, breaks = brks, main = "GAM")
plot(rf_proj_fav, col = clrs, breaks = brks, main = "RF")
plot(gbm_proj_fav, col = clrs, breaks = brks, main = "GBM")

par(mfrow = c(1, 3))
plot(maxent_proj, col = clrs, breaks = brks, main = "Maxent")
plot(glm_proj_fav, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_fav, col = clrs, breaks = brks, main = "GAM")

# ANALYSE ENVIRONMENTAL DISSIMILARITY ####

?mess  # multivariate environmental similarity surface
mess_proj <- mess(x = layers_proj, v = na.omit(getValues(layers_mod)))
mess_proj
plot(mess_proj, col = clrs, main = "MESS")  # negative values indicate environmental dissimilarity from the reference region
# in those places, our model projections are less reliable


# MODEL ENSEMBLES ####

# we can weigh each model according to its cross-validation result

# import the cross-validation results saved at the end of the previous practical:
crossval <- read.csv(paste0("./output/niche_model/crossval_", str_remove(myspecies, "Linum "), ".csv"))
crossval

# compute mean performance across the cross-validation folds:
crossval_means <- sapply(crossval, mean, na.rm = TRUE)
crossval_means

# choose e.g. two metrics (one discrimination and one calibration metric) to select models for the ensemble:
crossval_AUC <- crossval_means[grep("AUC", names(crossval_means))]
crossval_AUC
crossval_MCS <- crossval_means[grep("MCS", names(crossval_means))]
crossval_MCS
names(crossval_AUC) <- sapply(strsplit(names(crossval_AUC), "_"), `[`, 1)
names(crossval_MCS) <- sapply(strsplit(names(crossval_MCS), "_"), `[`, 1)
crossval_AUC
crossval_MCS

# set the acceptable performance thresholds:
thresh_AUC <- 0.7
thresh_MCS <- 0.5

# select models based on these performance thresholds & visual inspection of initial models (remove random forest):
mods_sel_AUC <- crossval_AUC[which(crossval_AUC >= thresh_AUC)]
mods_sel_MCS <- crossval_MCS[which(crossval_MCS >= thresh_MCS & crossval_MCS <= 1 + thresh_MCS)]
mods_sel_AUC
mods_sel_MCS

mods_sel_AUC <- mods_sel_AUC[-6] #remove random forest
mods_sel_MCS <- mods_sel_MCS[-6] #remove random forest


# make a raster stack of the prediction maps:
pred_stack <- stack(glm_fav, gam_fav, maxent_pred, maxnet_pred, gbm_fav)
names(pred_stack) <- c("glm", "gam", "maxent", "maxnet", "gbm")
plot(pred_stack, col = clrs, breaks = brks)


# conditional average (and variance): include only models that pass AUC and MCS thresholds

mods_sel_stack <- subset(pred_stack, intersect(names(mods_sel_AUC), names(mods_sel_MCS)))
plot(mods_sel_stack, col = clrs, breaks = brks)


par(mfrow = c(2, 1), mar = c(2, 2, 2, 1))

# average predictions of selected models:
pred_cond_avg <- mean(mods_sel_stack)
pred_cond_avg
plot(pred_cond_avg, col = clrs, breaks = brks, main = "Prediction mean")

# variance of the predictions of selected models:
pred_cond_avg_var <- calc(mods_sel_stack, var, use = "pairwise.complete.obs")
pred_cond_avg_var
plot(pred_cond_avg_var, col = clrs, main = "Prediction variance")
# prediction mean is more reliable where prediction variance is lower


# weighted average (and variance): weigh predictions by AUC (or another metric)

names(mods_sel_AUC)
names(mods_sel_stack)  # check that the names are in the same order!
pred_weight_avg <- weighted.mean(mods_sel_stack, w = mods_sel_AUC)
pred_cond_avg_var <- calc(mods_sel_stack, var, use = "pairwise.complete.obs")
plot(pred_weight_avg, col = clrs, breaks = brks, main = "AUC-weighted mean")
plot(pred_cond_avg_var, col = clrs, main = "Prediction variance")


# PROJECT TO DIFFERENT TIME PERIODS ####

# you can download future variables with 'sdmpredictors' package
# and from many online sources, see presentation "04_predictor_variables.pdf"
# let's take a peek at what's available in 'sdmpredictors' right now:

## FUTURE

##future_layers <- list_layers_future()
##head(future_layers)
##names(future_layers)
##unique(future_layers$dataset_code)  # currently only 2 datasets have projected future values

##future_worldclim <- subset(future_layers, dataset_code == "WorldClim")

# see which layers/scenarios are available for future WorldClim:
##unique(future_worldclim[ , c("model", "scenario", "year", "version")])


# select climate model(s), emissions scenario(s) and year(s) to project the models:
##unique(future_worldclim[ , c("model", "scenario", "year")])

##unique(subset(future_worldclim, model == "CCSM4" & scenario == "rcp26" & year == 2050))

##options(sdmpredictors_datadir = "./data/climate/future")

# example:
##vars_sel_2050_rcp26 <- subset(future_worldclim, model == "CCSM4" & scenario == "rcp26" & year == 2050)[1:19, ]$layer_code

##layers_2050_rcp26 <- load_layers(layercodes = vars_sel_2050_rcp26, rasterstack = FALSE)

##layers_f <- stack(layers_2050_rcp26)
##plot(layers)

##layers_2050_rcp26 <- raster::crop(layers_f, proj_extent)
##plot(layers_2050_rcp26, col = clrs)

##names(layers_2050_rcp26) <- str_replace(names(layers_2050_rcp26), "cc26_2050", "lonlat")

##layers_2050_rcp26_f <- stack(layers_2050_rcp26, layers_proj[[1]])
##plot(layers_2050_rcp26_f, col = clrs)

# then you can use these layers to project your models, as we did before when projecting to a different region

# predict with each model to the projection layers:
##bioclim_proj_f <- predict(layers_2050_rcp26_f, bioclim_mod)
##domain_proj_f <- predict(layers_2050_rcp26_f[[names(domain_mod@presence)]], domain_mod)  # takes time!
##maxent_proj_f <- predict(layers_2050_rcp26_f, maxent_mod, type = "cloglog")
##maxnet_proj_f <- predict(layers_2050_rcp26_f, maxnet_mod, type = "cloglog")
##glm_proj_f <- predict(layers_2050_rcp26_f, glm_mod, type = "response")
##gam_proj_f <- predict(layers_2050_rcp26_f, gam_mod, type = "response")
##rf_proj_f <- 1 - predict(layers_2050_rcp26_f, rf_mod, type = "prob")
##gbm_proj_f <- predict(layers_2050_rcp26_f, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
##glm_proj_fav_f <- Fav(pred = glm_proj_f, sample.preval = prevalence(glm_mod$y))
##gam_proj_fav_f <- Fav(pred = gam_proj_f, sample.preval = prevalence(gam_mod$y))
##rf_proj_fav_f <- Fav(pred = rf_proj_f, sample.preval = prevalence(rf_mod$y))
##gbm_proj_fav_f <- Fav(pred = gbm_proj_f, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):
##par(mfrow = c(4, 2))

##plot(bioclim_proj_f, col = clrs, breaks = brks, main = "Bioclim")
##plot(domain_proj_f, col = clrs, breaks = brks, main = "Domain")
##plot(maxent_proj_f, col = clrs, breaks = brks, main = "Maxent")
##plot(maxnet_proj_f, col = clrs, breaks = brks, main = "Maxnet")
##plot(glm_proj_fav_f, col = clrs, breaks = brks, main = "GLM")
##plot(gam_proj_fav_f, col = clrs, breaks = brks, main = "GAM")
##plot(rf_proj_fav_f, col = clrs, breaks = brks, main = "RF")
##plot(gbm_proj_fav_f, col = clrs, breaks = brks, main = "GBM")

## PAST

past_layers <- list_layers_paleo()
head(past_layers)
names(past_layers)
unique(past_layers$dataset_code)  # currently only 2 datasets have projected future values

past_worldclim <- subset(past_layers, dataset_code == "WorldClim")

# see which layers/scenarios are available for future WorldClim:
unique(past_worldclim[ , c("model_name", "epoch", "version")])
#LGM = 21000 years before present
#Mid Holocene = 6000 years before present

# select climate model(s), emissions scenario(s) and year(s) to project the models:
unique(subset(past_worldclim, model_name == "lgm CCSM4" | model_name == "holo CCSM4"))

options(sdmpredictors_datadir = "./data/niche_model/climate/past")

# example:
vars_sel_past <- subset(past_worldclim, model_name == "lgm CCSM4" | model_name == "holo CCSM4")[1:38, ]$layer_code

layers_past <- load_layers(layercodes = vars_sel_past, rasterstack = FALSE)

layers_past <- stack(layers_past)

layers_lgm <- raster::subset(layers_past, grep('cclgm', names(layers_past), value = T))
layers_hol <- raster::subset(layers_past, grep('ccmid', names(layers_past), value = T))

layers_lgm <- raster::crop(layers_lgm, proj_extent)
layers_lgm$WC_alt_lonlat <- layers_proj$WC_alt_lonlat #add altitude
plot(layers_lgm, col = clrs)

layers_hol <- raster::crop(layers_hol, proj_extent)
layers_hol$WC_alt_lonlat <- layers_proj$WC_alt_lonlat #add altitude
plot(layers_hol, col = clrs)

names(layers_lgm) <- str_replace(names(layers_lgm), "cclgm", "lonlat")
names(layers_hol) <- str_replace(names(layers_hol), "ccmid", "lonlat")


## Then also add LIG (Last Interglacial: 140-120000 years before present) directly downloaded from worldclim

files.lig <- list.files(path = "./data/niche_model/climate/past/lig_30s_bio",  pattern = "\\.bil$", full.names = TRUE)

r_name <-  list.files(path = "./data/niche_model/climate/past/lig_30s_bio",  pattern = "\\.bil$", full.names = FALSE) %>%
  str_replace("lig_30s_", "WC_") %>%
  str_replace("bio_", "bio") %>%
  str_replace(".bil", "_lonlat")

raster.list <- list() # to save raster values

for(i in 1:length(files.lig)){
  temp <- raster(files.lig[i])
  names(temp) <- r_name[i]
  raster.list[[i]] <- temp
}

layers_lig <- stack(raster.list)
layers_lig <- raster::crop(layers_lig, proj_extent)

alt.rsmp <- raster::crop(layers_proj$WC_alt_lonlat, proj_extent)
alt.rsmp <- resample(alt.rsmp, layers_lig[[1]]) #make sure raster layers for land surface and climate have same resolution

layers_lig$WC_alt_lonlat <- alt.rsmp
plot(layers_lig, col = clrs)


## Then also add PLIOCENE (3Mya years before present) directly downloaded from https://www.ecoclimate.org/downloads/

### Pliocene layers need to be cropped also for land surface so I use bio1 layer for Pliocene retrieved at http://www.paleoclim.org/
### Paleoclim would have been a better dataset to use directly but it did not have some layers I needed (e.g. bio2 and 3)

## Land
land.plio <- raster("./data/niche_model/climate/past/pliocene/M2_v1_r2_5m/bio_1.tif")
land.plio <- raster::crop(land.plio, proj_extent)

plot(land.plio)

## Climate
files.plio <- list.files(path = "./data/niche_model/climate/past/pliocene",  pattern = "\\.bil$", full.names = TRUE)

r_name <-  list.files(path = "./data/niche_model/climate/past/pliocene",  pattern = "\\.bil$", full.names = FALSE) %>%
  str_replace("bio..baseline_Historical.1900.1949...CCSM_Plio.3Ma._", "WC_") %>%
  str_replace(".bil", "_lonlat")

raster.list <- list() # to save raster values

for(i in 1:length(files.plio)){
  temp <- raster::crop(raster(files.plio[i]), proj_extent)
  land.plio1 <- resample(land.plio, temp) #make sure raster layers for land surface and climate have same resolution
  temp <- mask(temp, land.plio1)
  names(temp) <- r_name[i]
  raster.list[[i]] <- temp
}

layers_plio <- stack(raster.list)

## Altitude
alt.rsmp <- raster::crop(layers_proj$WC_alt_lonlat, proj_extent)
alt.rsmp <- resample(alt.rsmp, layers_plio[[1]]) #make sure raster layers for land surface and climate have same resolution

layers_plio$WC_alt_lonlat <- alt.rsmp #add altitude
plot(layers_plio, col = clrs)


### PLIOCENE 3.3 Mya
# predict with each model to the projection layers:
maxent_proj_plio <- predict(layers_plio, maxent_mod)
maxnet_proj_plio <- predict(layers_plio, maxnet_mod, type = "cloglog")
glm_proj_plio <- predict(layers_plio, glm_mod, type = "response")
gam_proj_plio <- predict(layers_plio, gam_mod, type = "response")
gbm_proj_plio <- predict(layers_plio, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_plioav_plio <- Fav(pred = glm_proj_plio, sample.preval = prevalence(glm_mod$y))
gam_proj_plioav_plio <- Fav(pred = gam_proj_plio, sample.preval = prevalence(gam_mod$y))
gbm_proj_plioav_plio <- Fav(pred = gbm_proj_plio, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(3, 2))

plot(maxent_proj_plio, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj_plio, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_plioav_plio, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_plioav_plio, col = clrs, breaks = brks, main = "GAM")
plot(gbm_proj_plioav_plio, col = clrs, breaks = brks, main = "GBM")


### LIG
# predict with each model to the projection layers:
maxent_proj_lig <- predict(layers_lig, maxent_mod, type = "cloglog")
maxnet_proj_lig <- predict(layers_lig, maxnet_mod, type = "cloglog")
glm_proj_lig <- predict(layers_lig, glm_mod, type = "response")
gam_proj_lig <- predict(layers_lig, gam_mod, type = "response")
gbm_proj_lig <- predict(layers_lig, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_ligav_lig <- Fav(pred = glm_proj_lig, sample.preval = prevalence(glm_mod$y))
gam_proj_ligav_lig <- Fav(pred = gam_proj_lig, sample.preval = prevalence(gam_mod$y))
gbm_proj_ligav_lig <- Fav(pred = gbm_proj_lig, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(3, 2))

plot(maxent_proj_lig, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj_lig, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_ligav_lig, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_ligav_lig, col = clrs, breaks = brks, main = "GAM")
plot(gbm_proj_ligav_lig, col = clrs, breaks = brks, main = "GBM")


### LGM
# predict with each model to the projection layers:
maxent_proj_lgm <- predict(layers_lgm, maxent_mod, type = "cloglog")
maxnet_proj_lgm <- predict(layers_lgm, maxnet_mod, type = "cloglog")
glm_proj_lgm <- predict(layers_lgm, glm_mod, type = "response")
gam_proj_lgm <- predict(layers_lgm, gam_mod, type = "response")
gbm_proj_lgm <- predict(layers_lgm, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_lgmav_lgm <- Fav(pred = glm_proj_lgm, sample.preval = prevalence(glm_mod$y))
gam_proj_lgmav_lgm <- Fav(pred = gam_proj_lgm, sample.preval = prevalence(gam_mod$y))
gbm_proj_lgmav_lgm <- Fav(pred = gbm_proj_lgm, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(3, 2))

plot(maxent_proj_lgm, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj_lgm, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_lgmav_lgm, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_lgmav_lgm, col = clrs, breaks = brks, main = "GAM")
plot(gbm_proj_lgmav_lgm, col = clrs, breaks = brks, main = "GBM")

## HOLOCENE
# predict with each model to the projection layers:
maxent_proj_hol <- predict(layers_hol, maxent_mod, type = "cloglog")
maxnet_proj_hol <- predict(layers_hol, maxnet_mod, type = "cloglog")
glm_proj_hol <- predict(layers_hol, glm_mod, type = "response")
gam_proj_hol <- predict(layers_hol, gam_mod, type = "response")
gbm_proj_hol <- predict(layers_hol, gbm_mod, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_holav_hol <- Fav(pred = glm_proj_hol, sample.preval = prevalence(glm_mod$y))
gam_proj_holav_hol <- Fav(pred = gam_proj_hol, sample.preval = prevalence(gam_mod$y))
gbm_proj_holav_hol <- Fav(pred = gbm_proj_hol, sample.preval = prevalence(gbm_mod$data$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(3, 2))

plot(maxent_proj_hol, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj_hol, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_holav_hol, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_holav_hol, col = clrs, breaks = brks, main = "GAM")
plot(gbm_proj_holav_hol, col = clrs, breaks = brks, main = "GBM")



# Mask prediction rasters for MESS < -50

mess_proj.mask <- mess_proj
mess_proj.mask[mess_proj.mask < -30] <- NA

maxent_proj.mskd <- mask(maxent_proj, mess_proj.mask)

mess_proj.mask.rsmp <- resample(mess_proj.mask, maxent_proj_plio)
maxent_proj_plio.mskd <- mask(maxent_proj_plio, mess_proj.mask.rsmp)

mess_proj.mask.rsmp <- resample(mess_proj.mask, maxent_proj_lig)
maxent_proj_lig.mskd <- mask(maxent_proj_lig, mess_proj.mask.rsmp)

maxent_proj_lgm.mskd <- mask(maxent_proj_lgm, mess_proj.mask)

maxent_proj_hol.mskd <- mask(maxent_proj_hol, mess_proj.mask)




world_map <- map_data("world")

eu_map <- ggplot(legend = FALSE) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray") +
  coord_equal(xlim = c(-20, 60), ylim = c(25, 60)) +
  theme_minimal()


maxent_proj.df <- as.data.frame(as(maxent_proj.mskd, "SpatialPixelsDataFrame"))

mapp <- eu_map +
  geom_tile(data = maxent_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "black", high = "white") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))


maxent_proj_hol.df <- as.data.frame(as(maxent_proj_hol.mskd, "SpatialPixelsDataFrame"))

maph <- eu_map +  
  geom_tile(data = maxent_proj_hol.df, aes(x = x, y = y, fill = layer), alpha = 0.8) +  
  borders("world", col = "black") +
  scale_fill_gradient(low = "black", high = "white") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))


maxent_proj_lgm.df <- as.data.frame(as(maxent_proj_lgm.mskd, "SpatialPixelsDataFrame"))

maplg <- eu_map +  
  geom_tile(data = maxent_proj_lgm.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") + 
  scale_fill_gradient(low = "black", high = "white") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))


maxent_proj_lig.df <- as.data.frame(as(maxent_proj_lig.mskd, "SpatialPixelsDataFrame"))

mapli <- eu_map +  
  geom_tile(data = maxent_proj_lig.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") + 
  scale_fill_gradient(low = "black", high = "white") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))


maxent_proj_plio.df <- as.data.frame(as(maxent_proj_plio.mskd, "SpatialPixelsDataFrame"))

mappl <- eu_map +  
  geom_tile(data = maxent_proj_plio.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") + 
  scale_fill_gradient(low = "black", high = "white") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))


mess_proj.df <- as.data.frame(as(mess_proj, "SpatialPixelsDataFrame"))
mess_proj.df$mess <- if_else(mess_proj.df$mess == Inf, NaN, mess_proj.df$mess)

mapms <- eu_map +
  geom_tile(data = mess_proj.df, aes(x = x, y = y, fill = mess)) + 
  borders("world", col = "black") +
  scale_fill_gradient2(low = "black", mid = "grey90", high = "white", midpoint = -30, na.value = "transparent") +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))




library(patchwork)

mapp + ((maph + maplg) / (mapli + mappl))


mapp
maph
maplg
mapli
mappl
mapms





















