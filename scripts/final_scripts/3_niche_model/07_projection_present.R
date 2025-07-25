# LIBS
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
library(sf)


# IMPORT DATA AND DEFINE SOME SETTINGS ####

# define the target species:
spp <- c("linum bienne")

# import models saved in previous practical:
load(paste0("./data/niche_model/models_", "lblus", ".RData"))

# define colours and colour breaks for the maps:
clrs <- hcl.colors(10)
brks <- seq(0, 1, by = 0.1)


# PLOT VARIABLE RESPONSE CURVES ####

# before you project your models, you may want to take a look at the response curves, to see if they make ecological sense
# but mind that individual response curves may be affected by interactions with other variables

response(maxent_mod.lb)
plot(maxnet_mod.lb)
plotmo(glm_mod.lb, trace = 1, all1 = TRUE)  # see the "plots" pane
plotmo(gam_mod.lb, trace = 1, all1 = TRUE)


# PROJECT THE MODELS TO A DIFFERENT REGION ####

# import the occurrences data:
sppd <- read.csv("./data/niche_model/species_occurrence/linumspp_1700today_geoclean.csv")

sppd <- sppd %>% filter(species %in% spp) %>% dplyr::select(-c(X.1, X)) %>% as.data.frame()

sppd.lb <- sppd %>% filter(species %in% spp[1]) %>% as.data.frame()

# convert to a spatial object:
occurrences_spatial.lb <- sppd.lb

coordinates(occurrences_spatial.lb) <- sppd.lb[ , c("decimalLongitude", "decimalLatitude")]
crs(occurrences_spatial.lb) <- "+proj=longlat"
par(mfrow = c(1, 1))
plot(occurrences_spatial.lb)

# import the countries map:
countries <- st_read("./data/niche_model/countries/world_countries.shp") 

plot(occurrences_spatial.lb, col = "blue")
plot(sf::as_Spatial(countries), 
     xlim = range(sppd.lb$decimalLongitude), ylim = range(sppd.lb$decimalLatitude),
     border = "black", col = "transparent", 
     lwd = 2, add = TRUE)
map.axes()


# import the complete (global coverage) layers:
layers <- raster::stack(list.files(path = "./data/niche_model/climate/wc2.1_10m_bio/", pattern = "*_crpd.tif$", full.names = TRUE))
#plot(layers)

# now choose the projection extent, i.e. the region where you want to project your model predictions
# you can use e.g. the wider area that contains species occurrence points (including those that were left out of the modelling region):
proj_extent <- extent(c(-24, 79, 22, 70))
# or you can choose a specific spatial window that interests you: #DIDN't do it, but might use a buffer or draw a polygon for my species from Diederichsen paper
#proj_extent <- extent(c(-15, 10, 30, 45))

# crop the global layers to the projection extent:
layers_proj <- crop(layers, proj_extent)
layers_proj <- raster::subset(layers_proj, vars_sel.lb)
#plot(layers_proj, col = clrs)

# predict with each model to the projection layers:
maxent_proj.lb <- predict(layers_proj, maxent_mod.lb, type = "cloglog")
maxnet_proj.lb <- predict(layers_proj, maxnet_mod.lb, type = "cloglog")
glm_proj.lb <- predict(layers_proj, glm_mod.lb, type = "response")
gam_proj.lb <- predict(layers_proj, gam_mod.lb, type = "response")

# convert probability predictions (from presence/absence models) to favourability:
glm_proj_fav.lb <- Fav(pred = glm_proj.lb, sample.preval = prevalence(glm_mod.lb$y))
gam_proj_fav.lb <- Fav(pred = gam_proj.lb, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(2,2))
plot(maxent_proj.lb, col = clrs, breaks = brks, main = "Maxent")
plot(maxnet_proj.lb, col = clrs, breaks = brks, main = "Maxnet")
plot(glm_proj_fav.lb, col = clrs, breaks = brks, main = "GLM")
plot(gam_proj_fav.lb, col = clrs, breaks = brks, main = "GAM")


# ANALYSE ENVIRONMENTAL DISSIMILARITY ####

# negative values indicate environmental dissimilarity from the reference region
# in those places, our model projections are less reliable

dev.off()

?mess  # multivariate environmental similarity surface
mess_proj.lb <- mess(x = layers_proj, v = na.omit(getValues(layers_mod.lb[[vars_sel.lb]])))
plot(mess_proj.lb, col = clrs, main = "MESS")  


# MODEL ENSEMBLES ####

# make a raster stack of the prediction maps:
pred_stack.lb <- stack(glm_proj_fav.lb, gam_proj_fav.lb, maxnet_proj.lb, maxent_proj.lb)
names(pred_stack.lb) <- c("glm", "gam", "maxnet", "maxent")
plot(pred_stack.lb, col = clrs, breaks = brks)

# conditional average (and variance): include only models that pass AUC and MCS thresholds
par(mfrow = c(2, 1), mar = c(2, 2, 2, 1))

# average predictions of selected models:
## LB
pred_cond_avg.lb <- mean(pred_stack.lb)
pred_cond_avg.lb
plot(pred_cond_avg.lb, col = clrs, breaks = brks, main = "Prediction mean")

# variance of the predictions of selected models:
pred_cond_avg_var.lb <- calc(pred_stack.lb, var, use = "pairwise.complete.obs")
pred_cond_avg_var.lb
plot(pred_cond_avg_var.lb, col = clrs, main = "Prediction variance")
# prediction mean is more reliable where prediction variance is lower

###### THRESHOLDING FOR L. BIENNE ######
msk <- mess_proj.lb
values(msk)[values(msk) < -20] <- NA
values(msk)[values(msk) == Inf] <- NA

writeRaster(mess_proj.lb, "./data/niche_model/MESS.tif", overwrite = TRUE)
writeRaster(msk, "./data/niche_model/MESS_mask.tif", overwrite = TRUE)


## Maxent
maxent_proj.lb.thr <- maxent_proj.lb
values(maxent_proj.lb.thr)[values(maxent_proj.lb.thr) < 0.3] <- NA
maxent_proj.lb.thr1 <- raster::mask(maxent_proj.lb.thr, msk) #crop areas of high env. dissimilarity

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
plot(maxent_proj.lb)
plot(maxent_proj.lb.thr1)

maxent_proj.lb.thr1 <- maxent_proj.lb.thr1 < 0.3 #distribution raster L. bienne
values(maxent_proj.lb.thr1)[values(maxent_proj.lb.thr1) == 0] <- 1
plot(maxent_proj.lb.thr1)

writeRaster(maxent_proj.lb, "./data/niche_model/maxent_preds_lb.tif", overwrite = TRUE)
writeRaster(maxent_proj.lb.thr1, "./data/niche_model/maxent_preds_lb_thr.tif", overwrite = TRUE)

## GAM
gam_proj.lb.thr <- gam_proj_fav.lb
values(gam_proj.lb.thr)[values(gam_proj.lb.thr) < 0.3] <- NA
gam_proj.lb.thr1 <- raster::mask(gam_proj.lb.thr, msk) #crop areas of high env. dissimilarity

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
plot(gam_proj_fav.lb)
plot(gam_proj.lb.thr1)

gam_proj.lb.thr1 <- gam_proj.lb.thr1 < 0.3 #distribution raster L. bienne
values(gam_proj.lb.thr1)[values(gam_proj.lb.thr1) == 0] <- 1
plot(gam_proj.lb.thr1)

writeRaster(gam_proj_fav.lb, "./data/niche_model/gam_preds_lb.tif", overwrite = TRUE)
writeRaster(gam_proj.lb.thr1, "./data/niche_model/gam_preds_lb_thr.tif", overwrite = TRUE)


save.image(paste0("./data/niche_model/models_", "lblus", ".RData"))
