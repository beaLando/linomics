# LIBS
library(dismo)
library(plotmo)
library(maxnet)
library(gam)
library(randomForest)
library(gbm)
library(raster)
library(pastclim)
library(maps)
library(fuzzySim)
library(tidyverse)

# IMPORT DATA AND DEFINE SOME SETTINGS ####

# define the target species:
spp <- c("linum bienne")

# import models saved in previous steps:
load(paste0("./data/niche_model/models_", "lblus", ".RData"))

# define colours and colour breaks for the maps:
clrs <- hcl.colors(10)
brks <- seq(0, 1, by = 0.1)

proj_extent <- extent(c(-24, 79, 22, 70))
proj_extenT <- terra::ext(c(-24, 79, 22, 70))

################ PAST PROJECTIONS ##################
# Use pastclim package to download 
#Mid Holocene = 8k years before present (mh) > temperate forest in north and mediterranean forest in south of italy
#LGM = ca. 21K years years before present (lgm)
#In between LGM and LIG = ca. 70K years before present (alig)
#LIG = ca. 130K years before present (lig)
#Mid Pleistocene = ca. 700K years bp (mple) > extinction of last subtropical taxa in italian peninsula
#Early Pleistocene = ca. 1.3M years bp (eple) > subtropical taxa still present, but expansion of herbaceous species in italian peninsula
#Late Pliocene = ca. 2.7M years bp (lpli) > onset of more seasonal climate in Med. basin
#Early Pliocene = ca.4.5M years bp (epli) > subtropical climate

set_data_path(path_to_nc = "./data/niche_model/climate/past")

# help("Barreto2023")#only first time
# get_vars_for_dataset(dataset = "Barreto2023") #ok: bio4, bio11, and bio12 are included
# 
# ?download_dataset
# download_dataset(dataset = "Barreto2023", bio_variables = c("bio04", "bio11", "bio12"))
get_time_bp_steps(dataset = "Barreto2023")

midh <- region_slice(
  time_bp = -8000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(midh)

lgm <- region_slice(
  time_bp = -22000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(lgm)

alig <- region_slice(
  time_bp = -70000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(alig)

lig <- region_slice(
  time_bp = -130000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(lig)

mple <- region_slice(
  time_bp = -700000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(mple)

eple <- region_slice(
  time_bp = -1300000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(eple)

lpli <- region_slice(
  time_bp = -2700000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(lpli)

epli <- region_slice(
  time_bp = -4500000,
  ext = proj_extenT,
  bio_variables = c("bio04", "bio11", "bio12"),
  dataset = "Barreto2023")

plot(epli)

## PREDICTIONS
### EARLY PLIOCENE
vars_sel.lb
varnames(epli) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(epli) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
epli <- stack(epli)

# predict with each model to the projection layers:
maxent_proj_epli <- predict(epli, maxent_mod.lb, na.rm = TRUE)
gam_proj_epli <- predict(epli, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_epliav_epli <- Fav(pred = gam_proj_epli, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_epli, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_epliav_epli, col = clrs, breaks = brks, main = "GAM")


### LATE PLIOCENE
varnames(lpli) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(lpli) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
lpli <- stack(lpli)

# predict with each model to the projection layers:
maxent_proj_lpli <- predict(lpli, maxent_mod.lb, na.rm = TRUE)
gam_proj_lpli <- predict(lpli, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_lpliav_lpli <- Fav(pred = gam_proj_lpli, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_lpli, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_lpliav_lpli, col = clrs, breaks = brks, main = "GAM")


### EARLY PLEISTOCENE
varnames(eple) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(eple) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
eple <- stack(eple)

# predict with each model to the projection layers:
maxent_proj_eple <- predict(eple, maxent_mod.lb, na.rm = TRUE)
gam_proj_eple <- predict(eple, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_epleav_eple <- Fav(pred = gam_proj_eple, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_eple, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_epleav_eple, col = clrs, breaks = brks, main = "GAM")


### MID PLEISTOCENE
varnames(mple) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(mple) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
mple <- stack(mple)

# predict with each model to the projection layers:
maxent_proj_mple <- predict(mple, maxent_mod.lb, na.rm = TRUE)
gam_proj_mple <- predict(mple, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_mpleav_mple <- Fav(pred = gam_proj_mple, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_mple, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_mpleav_mple, col = clrs, breaks = brks, main = "GAM")


### LIG
varnames(lig) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(lig) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
lig <- stack(lig)

# predict with each model to the projection layers:
maxent_proj_lig <- predict(lig, maxent_mod.lb, na.rm = TRUE)
gam_proj_lig <- predict(lig, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_ligav_lig <- Fav(pred = gam_proj_lig, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_lig, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_ligav_lig, col = clrs, breaks = brks, main = "GAM")


### Between LIG and LGM
varnames(alig) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(alig) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
alig <- stack(alig)

# predict with each model to the projection layers:
maxent_proj_alig <- predict(alig, maxent_mod.lb, na.rm = TRUE)
gam_proj_alig <- predict(alig, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_aligav_alig <- Fav(pred = gam_proj_alig, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_alig, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_aligav_alig, col = clrs, breaks = brks, main = "GAM")


### LGM
varnames(lgm) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(lgm) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
lgm <- stack(lgm)

# predict with each model to the projection layers:
maxent_proj_lgm <- predict(lgm, maxent_mod.lb, na.rm = TRUE)
gam_proj_lgm <- predict(lgm, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_lgmav_lgm <- Fav(pred = gam_proj_lgm, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_lgm, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_lgmav_lgm, col = clrs, breaks = brks, main = "GAM")


### MID-HOLOCENE
varnames(midh) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
names(midh) <- c("bio_4_crpd", "bio_11_crpd", "bio_12_crpd")
midh <- stack(midh)

# predict with each model to the projection layers:
maxent_proj_midh <- predict(midh, maxent_mod.lb, na.rm = TRUE)
gam_proj_midh <- predict(midh, gam_mod.lb, type = "response", na.rm = TRUE)
gam_proj_midhav_midh <- Fav(pred = gam_proj_midh, sample.preval = prevalence(gam_mod.lb$y))

# map the extrapolated model predictions (model projections):
par(mfrow = c(1, 2))

plot(maxent_proj_midh, col = clrs, breaks = brks, main = "Maxent")
plot(gam_proj_midhav_midh, col = clrs, breaks = brks, main = "GAM")



### ALL - mask projections with mess
mskr <- resample(msk, maxent_proj_midh)

## Maxent
maxent_proj.lb.mskrd <- raster::mask(maxent_proj.lb, msk) #crop areas of high env. dissimilarity
maxent_proj_epli.mskrd <- mask(maxent_proj_epli, mskr)
maxent_proj_lpli.mskrd <- mask(maxent_proj_lpli, mskr)
maxent_proj_eple.mskrd <- mask(maxent_proj_eple, mskr)
maxent_proj_mple.mskrd <- mask(maxent_proj_mple, mskr)
maxent_proj_lig.mskrd <- mask(maxent_proj_lig, mskr)
maxent_proj_alig.mskrd <- mask(maxent_proj_alig, mskr)
maxent_proj_lgm.mskrd <- mask(maxent_proj_lgm, mskr)
maxent_proj_hol.mskrd <- mask(maxent_proj_midh, mskr)

par(mfrow = c(3, 3))
plot(maxent_proj_epli.mskrd,zlim=c(0,1))
plot(maxent_proj_lpli.mskrd,zlim=c(0,1))
plot(maxent_proj_eple.mskrd,zlim=c(0,1))
plot(maxent_proj_mple.mskrd,zlim=c(0,1))
plot(maxent_proj_lig.mskrd,zlim=c(0,1))
plot(maxent_proj_alig.mskrd,zlim=c(0,1))
plot(maxent_proj_lgm.mskrd,zlim=c(0,1))
plot(maxent_proj_hol.mskrd,zlim=c(0,1))
plot(maxent_proj.lb.mskrd,zlim=c(0,1))

## GAM
gam_proj.lb.mskrd <- raster::mask(gam_proj.lb, msk) #crop areas of high env. dissimilarity
gam_proj_epli.mskrd <- mask(gam_proj_epli, mskr)
gam_proj_lpli.mskrd <- mask(gam_proj_lpli, mskr)
gam_proj_eple.mskrd <- mask(gam_proj_eple, mskr)
gam_proj_mple.mskrd <- mask(gam_proj_mple, mskr)
gam_proj_lig.mskrd <- mask(gam_proj_lig, mskr)
gam_proj_alig.mskrd <- mask(gam_proj_alig, mskr)
gam_proj_lgm.mskrd <- mask(gam_proj_lgm, mskr)
gam_proj_hol.mskrd <- mask(gam_proj_midh, mskr)

par(mfrow = c(3, 3))
plot(gam_proj_epli.mskrd,zlim=c(0,1))
plot(gam_proj_lpli.mskrd,zlim=c(0,1))
plot(gam_proj_eple.mskrd,zlim=c(0,1))
plot(gam_proj_mple.mskrd,zlim=c(0,1))
plot(gam_proj_lig.mskrd,zlim=c(0,1))
plot(gam_proj_alig.mskrd,zlim=c(0,1))
plot(gam_proj_lgm.mskrd,zlim=c(0,1))
plot(gam_proj_hol.mskrd,zlim=c(0,1))
plot(gam_proj.lb.mskrd,zlim=c(0,1))

# SAVE as RASTERS
writeRaster(maxent_proj_epli, "./data/niche_model/maxent_preds_lb_1_earlyPlio.tif", overwrite = TRUE)
writeRaster(maxent_proj_lpli, "./data/niche_model/maxent_preds_lb_2_latePlio.tif", overwrite = TRUE)
writeRaster(maxent_proj_eple, "./data/niche_model/maxent_preds_lb_3_earlyPlei.tif", overwrite = TRUE)
writeRaster(maxent_proj_mple, "./data/niche_model/maxent_preds_lb_4_midPlei.tif", overwrite = TRUE)
writeRaster(maxent_proj_lig, "./data/niche_model/maxent_preds_lb_5_lig.tif", overwrite = TRUE)
writeRaster(maxent_proj_alig, "./data/niche_model/maxent_preds_lb_6_lig-lgm.tif", overwrite = TRUE)
writeRaster(maxent_proj_lgm, "./data/niche_model/maxent_preds_lb_7_lgm.tif", overwrite = TRUE)
writeRaster(maxent_proj_midh, "./data/niche_model/maxent_preds_lb_8_holo.tif", overwrite = TRUE)

writeRaster(gam_proj_epli, "./data/niche_model/gam_preds_lb_1_earlyPlio.tif", overwrite = TRUE)
writeRaster(gam_proj_lpli, "./data/niche_model/gam_preds_lb_2_latePlio.tif", overwrite = TRUE)
writeRaster(gam_proj_eple, "./data/niche_model/gam_preds_lb_3_earlyPlei.tif", overwrite = TRUE)
writeRaster(gam_proj_mple, "./data/niche_model/gam_preds_lb_4_midPlei.tif", overwrite = TRUE)
writeRaster(gam_proj_lig, "./data/niche_model/gam_preds_lb_5_lig.tif", overwrite = TRUE)
writeRaster(gam_proj_alig, "./data/niche_model/gam_preds_lb_6_lig-lgm.tif", overwrite = TRUE)
writeRaster(gam_proj_lgm, "./data/niche_model/gam_preds_lb_7_lgm.tif", overwrite = TRUE)
writeRaster(gam_proj_midh, "./data/niche_model/gam_preds_lb_8_holo.tif", overwrite = TRUE)


# PLOTTTTS as MAPPPS
# ## Read in rasters of predictions
# maxent_proj.lb <- raster("./data/niche_model/maxent_preds_lb.tif")

## Plot each epoch
world_map <- map_data("world")

eu_map <- ggplot(legend = FALSE) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "white") +
  coord_equal(xlim = c(-20, 60), ylim = c(25, 60)) +
  theme_minimal()

### Present
gam_proj.df <- as.data.frame(as(gam_proj.lb.mskrd, "SpatialPixelsDataFrame"))

mapp <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### Holocene
gam_proj.df <- as.data.frame(as(gam_proj_hol.mskrd, "SpatialPixelsDataFrame"))

maph <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### LGM
gam_proj.df <- as.data.frame(as(gam_proj_lgm.mskrd, "SpatialPixelsDataFrame"))

maplgm <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### LGM-LIG
gam_proj.df <- as.data.frame(as(gam_proj_alig.mskrd, "SpatialPixelsDataFrame"))

mapali <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### LIG
gam_proj.df <- as.data.frame(as(gam_proj_lig.mskrd, "SpatialPixelsDataFrame"))

maplig <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### Mid-PLEISTOCENE
gam_proj.df <- as.data.frame(as(gam_proj_mple.mskrd, "SpatialPixelsDataFrame"))

mapmple <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### EARLY-PLEISTOCENE
gam_proj.df <- as.data.frame(as(gam_proj_eple.mskrd, "SpatialPixelsDataFrame"))

mapeple <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### LATE PLIOCENE
gam_proj.df <- as.data.frame(as(gam_proj_lpli.mskrd, "SpatialPixelsDataFrame"))

maplpli <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### EARLY PLIOCENE
gam_proj.df <- as.data.frame(as(gam_proj_epli.mskrd, "SpatialPixelsDataFrame"))

mapepli <- eu_map +
  geom_tile(data = gam_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  borders("world", col = "black") +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent", limits=c(0,1)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Longitude", y = "Latitude", fill = "Favourability") +
  theme(legend.position = c("right"))

### ALL TOGETHER
# library(patchwork)
# mapp + ((maph + maplgm + mapali + maplig) / (mapmple + mapeple + maplpli + mapepli)) + plot_layout(guides="collect")

mapp
maph
maplgm
mapali
maplig
mapmple
mapeple
maplpli
mapepli