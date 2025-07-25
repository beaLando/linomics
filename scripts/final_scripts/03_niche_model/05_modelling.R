# R scripts for the course "Species distribution and ecological niche modelling in R"
# A. Marcia Barbosa (https://modtools.wordpress.com)


# don't forget to DO THIS FIRST:
# Session -> Set Working Directory -> To Source File Location


# LOAD PACKAGES ####
library(raster)
library(fuzzySim)
library(modEvA)
library(dismo)
library(maxnet)
library(gam)
library(randomForest)
library(gbm)
library(corrplot)
library(tidyverse)


# IMPORT DATA AND DEFINE SOME SETTINGS ####

spp <- c("linum bienne")

# import the data saved during the previous practical:
dat.lb <- read.csv(paste0("./data/niche_model/dat_", spp[1], ".csv"))
head(dat.lb)

# or (depending on which dataset you want to use):
#dat <- read.csv(paste0("../outputs/dat_grid_", myspecies, ".csv"))
layers_mod.lb <- readRDS(paste0("./data/niche_model/climate_mod_", spp[1], ".rds"))
plot(layers_mod.lb)

# define a colour palette and colour breaks, so that all maps of model predictions show a comparable gradient:
clrs <- hcl.colors(10)
brks <- seq(0, 1, by = 0.1)


# just to check when we are actually modelling presence data only, let's make a version of the 'layers_mod' RasterStack with pixel values only at the locations of the presence points:
##LB
names(dat.lb)
presence_centroids.lb <- dat.lb[dat.lb$presence == 1, ]
coordinates(presence_centroids.lb) <- presence_centroids.lb[ , c("x", "y")]
crs(presence_centroids.lb) <- "+proj=longlat"
plot(presence_centroids.lb)

layers_mod_presonly.lb <- raster::mask(layers_mod.lb, presence_centroids.lb)
plot(layers_mod_presonly.lb, col = clrs)
plot(layers_mod_presonly.lb[[1]], col = clrs)


# SELECT VARIABLES ####

# some bioclimatic variables were found to correlate poorly among present and past/future scenarios (https://doi.org/10.1371/journal.pone.0129037, https://doi.org/10.1016/j.gloplacha.2013.04.005); for actual work, you may want to look at these papers and see if this affects your study region and time periods

names(dat.lb)
vars.lb <- names(dat.lb)[grep("ele|bio", names(dat.lb))]
vars.lb  # check if OK for your dataset!

# keep only variables also present in pliocene datasets
vars.plio <- c("bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19", "bio_4", "bio_8", "bio_9")
vars.plio <- paste0(vars.plio, "_crpd")
vars.diff <- setdiff(vars.lb, vars.plio)

vars.lb <- vars.lb[-c(which(vars.lb %in% vars.diff))]

corrplot(cor(dat.lb[which(dat.lb$presence==1),5:23]), type = "upper", method = "number")
corrplot(cor(dat.lb[which(dat.lb$presence==1),c(5,8,11,18,13,7)]), type = "upper", method = "number")

keep <- colnames(dat.lb)[c(5,8,11,18,13,7)]

keep %in% vars.lb #all included

vars.lb <- vars.lb[which(vars.lb %in% keep)]

# many methods are affected by having too many variables and/or high correlations among them
# so, let's select a subset of not highly correlated variables to use with all models (so that they are comparable)
# (note this is just one possible way of selecting variables!)

vars_sel.lb <- corSelect(data = dat.lb, sp.cols = "presence", 
                         var.cols = vars.lb, 
                         cor.thresh = 0.6, 
                         select = "AIC")
vars_sel.lb
vars_sel.lb <- vars_sel.lb$selected.vars
vars_sel.lb

plot(layers_mod.lb[[vars_sel.lb]])


# PRESENCE/BACKGROUND MODELLING: MAXENT ####

## MAXENT
?maxent  # 'maxent' function in 'dismo' package
maxent_mod.lb <- maxent(x = layers_mod.lb[[vars_sel.lb]], p = presence_centroids.lb)

# let's look at the contents of the resulting object:
str(maxent_mod.lb)
nrow(maxent_mod.lb@presence)
nrow(maxent_mod.lb@absence)  # it creates absences from the pixels without presence points! so it does use absence data
# compute and map maxent predictions:
maxent_pred.lb <- predict(layers_mod.lb, maxent_mod.lb)
plot(maxent_pred.lb, col = clrs, breaks = brks, main = "Maxent")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)

## MAXNET
maxnet_mod.lb <- maxnet(p = dat.lb$presence, data = dat.lb[ , vars_sel.lb], maxnet.formula(p = dat.lb$presence, data = dat.lb[ , vars_sel.lb]))

# compute and map maxnet predictions:
maxnet_pred.lb <- predict(layers_mod.lb, maxnet_mod.lb, type = "cloglog")
plot(maxnet_pred.lb, col = clrs, breaks = brks, main = "Maxnet")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)


# PRESENCE/ABSENCE MODELLING ####

## GLM
# let's first make the formula required by the 'glm' function:
glm_form.lb <- as.formula(paste("presence ~", paste(vars_sel.lb, collapse = " + ")))
glm_form.lb

glm_mod.lb <- glm(formula = glm_form.lb, family = binomial, data = dat.lb)
summary(glm_mod.lb, correlation = TRUE)

glm_pred.lb <- predict(layers_mod.lb, glm_mod.lb, type = "response")  # "response" transforms the linear predictor with the appropriate link function, yielding results as probability values
plot(glm_pred.lb, col = clrs, breaks = brks, main = "GLM")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability:
glm_fav.lb <- Fav(pred = glm_pred.lb, sample.preval = prevalence(glm_mod.lb$y))
plot(glm_fav.lb, col = clrs, breaks = brks, main = "GLM favourability")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)

## GAM
gam_form.lb <- as.formula(paste("presence ~", paste0("s(", vars_sel.lb, ")", collapse = "+")))  # GAM with smoothing splines ('s')
gam_form.lb

gam_mod.lb <- gam(formula = gam_form.lb, family = binomial, data = dat.lb)
summary(gam_mod.lb)

gam_pred.lb <- predict(layers_mod.lb, gam_mod.lb, type = "response")
plot(gam_pred.lb, col = clrs, breaks = brks, main = "GAM")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability:
gam_fav.lb <- Fav(pred = gam_pred.lb, sample.preval = prevalence(gam_mod.lb$y))
plot(gam_fav.lb, col = clrs, breaks = brks, main = "GAM favourability")
plot(presence_centroids.lb, pch = ".", col = "red", add = TRUE)

# plot the raster map predictions together:
par(mar = c(2, 2, 2, 1), mfrow = c(2, 2))
plot(maxent_pred.lb, col = clrs, main = "Maxent")
plot(maxnet_pred.lb, col = clrs, main = "Maxnet")
plot(glm_fav.lb, col = clrs, main = "GLM fav")
plot(gam_fav.lb, col = clrs, main = "GAM fav")


# APPLY PREDICTIONS TO THE DATA TABLE ####
# (variables must have exact same names as in the models):
dat.lb$maxent_pred <- predict(maxent_mod.lb, dat.lb)
dat.lb$maxnet_pred <- as.vector(predict(maxnet_mod.lb, dat.lb, type = "cloglog"))
dat.lb$glm_pred <- predict(glm_mod.lb, newdata = dat.lb, type = "response")
dat.lb$glm_fav <- Fav(pred = dat.lb$glm_pred, sample.preval = prevalence(glm_mod.lb$y))
dat.lb$gam_pred <- predict(gam_mod.lb, newdata = dat.lb, type = "response")
dat.lb$gam_fav <- Fav(pred = dat.lb$gam_pred, sample.preval = prevalence(gam_mod.lb$y))
head(dat.lb)

pred_cols.lb <- c("maxent_pred", "maxnet_pred", "glm_fav", "gam_fav")


# map predictions from the data table:

# first we'll need to convert 'dat' to a spatial object:
dat_spatial.lb <- dat.lb
names(dat_spatial.lb)
coordinates(dat_spatial.lb) <- dat_spatial.lb[ , c("x", "y")]

# map all prediction columns in one window:
spplot(dat_spatial.lb, zcol = pred_cols.lb, col.regions = clrs, cex = 0.1)


# COMPARE PREDICTIONS OF DIFFERENT MODELS ####

# pair-wise scatterplots of predictions:
pairs(dat.lb[ , pred_cols.lb], pch = 20, cex = 0.1)  # takes time!

# correlations among predictions:
pred_corrs.lb <- cor(dat.lb[ , pred_cols.lb])
pred_corrs.lb
min(pred_corrs.lb, na.rm = TRUE) #0.74

# visual correlation matrix:
par(mfrow = c(1, 1))
corrplot(pred_corrs.lb, method = "ellipse", type = "upper", addCoef.col = "wheat3", addCoefasPercent = TRUE)

# save all current objects to disk (YOU WILL NEED THEM LATER!):
save.image(paste0("./data/niche_model/models_", "lblus", ".RData"))
