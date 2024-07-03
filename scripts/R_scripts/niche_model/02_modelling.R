# R scripts for the course "Species distribution and ecological niche modelling in R"
# A. Marcia Barbosa (https://modtools.wordpress.com)


# don't forget to DO THIS FIRST:
# Session -> Set Working Directory -> To Source File Location


# LOAD PACKAGES ####
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

myspecies <- c("Linum bienne")

# import the data saved during the previous practical:
dat <- read.csv(paste0("./data/niche_model/4modelling/", str_remove(myspecies, "Linum "), ".csv"))
# or (depending on which dataset you want to use):
#dat <- read.csv(paste0("../outputs/dat_grid_", myspecies, ".csv"))
head(dat)
layers_mod <- readRDS(paste0("./data/niche_model/4modelling/", str_remove(myspecies, "Linum "), ".rds"))
plot(layers_mod)


# define a colour palette and colour breaks, so that all maps of model predictions show a comparable gradient:
clrs <- hcl.colors(10)
brks <- seq(0, 1, by = 0.1)


# just to check when we are actually modelling presence data only, let's make a version of the 'layers_mod' RasterStack with pixel values only at the locations of the presence points:

names(dat)
presence_centroids <- dat[dat$presence == 1, ]
coordinates(presence_centroids) <- presence_centroids[ , c("x", "y")]
crs(presence_centroids) <- "+proj=longlat"
plot(presence_centroids)

layers_mod_presonly <- mask(layers_mod, presence_centroids)
plot(layers_mod_presonly, col = clrs)
plot(layers_mod_presonly[[1]], col = clrs)


# SELECT VARIABLES ####

# some bioclimatic variables were found to correlate poorly among present and past/future scenarios (https://doi.org/10.1371/journal.pone.0129037, https://doi.org/10.1016/j.gloplacha.2013.04.005); for actual work, you may want to look at these papers and see if this affects your study region and time periods

names(dat)
vars <- names(dat)[grep("alt|bio", names(dat))]
vars  # check if OK for your dataset!

# if you want to exclude particular variables:
#exclude <- names(dat)[grep("bio3|bio14|bio15", names(dat))]
#exclude  # check if OK for your dataset!
#vars <- setdiff(vars, exclude)
#vars  # check if OK for your dataset!

# many methods are affected by having too many variables and/or high correlations among them
# so, let's select a subset of not highly correlated variables to use with all models (so that they are comparable)
# (note this is just one possible way of selecting variables!)
?corSelect
vars_sel <- corSelect(data = dat, sp.cols = "presence", var.cols = vars, cor.thresh = 0.8, select = "AIC")
vars_sel
vars_sel <- vars_sel$selected.vars
vars_sel

vars_sel2 <- corSelect(data = dat, sp.cols = "presence", var.cols = vars[-1], cor.thresh = 0.8, select = "AIC") #no altitude
vars_sel2
vars_sel2 <- vars_sel2$selected.vars
vars_sel2

# PRESENCE-ONLY MODELLING: BIOCLIM ####

## with altitude
bioclim_mod_po <- bioclim(x = layers_mod_presonly[[vars_sel]], p = presence_centroids)  # presence-only
bioclim_mod <- bioclim(x = layers_mod[[vars_sel]], p = presence_centroids)
all.equal(bioclim_mod, bioclim_mod_po)  # TRUE, so we confirm that we get the same model whether or not we provide pixels without presence records, showing that Bioclim is a true presence-only method

bioclim_mod
plot(bioclim_mod)  # the axes represent only the first 2 variables in the model, but the response includes all variables in the model

bioclim_pred <- predict(layers_mod, bioclim_mod)

par(mar = c(2, 2, 1, 1))
plot(bioclim_pred, col = clrs, breaks = brks, main = "Bioclim")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

## without altitude
bioclim_mod_na <- bioclim(x = layers_mod[[vars_sel2]], p = presence_centroids)
bioclim_mod_na
plot(bioclim_mod_na)

bioclim_pred_na <- predict(layers_mod, bioclim_mod_na)

par(mar = c(2, 2, 1, 1))
plot(bioclim_pred_na, col = clrs, breaks = brks, main = "Bioclim")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


# PRESENCE-ONLY MODELLING: DOMAIN ####

## with altitude
?domain

domain_mod_po <- domain(x = layers_mod_presonly[[vars_sel]], p = presence_centroids)  # presence-only
domain_mod <- domain(x = layers_mod[[vars_sel]], p = presence_centroids)
all.equal(domain_mod, domain_mod_po)  # TRUE, so Domain also provides the same model whether or not we give it the values of the pixels without presence, which means it is a true presence-only method

domain_mod

domain_pred <- predict(layers_mod[[vars_sel]], domain_mod)

## without altitude
domain_mod_na <- domain(x = layers_mod[[vars_sel2]], p = presence_centroids)
domain_mod_na

domain_pred_na <- predict(layers_mod[[vars_sel2]], domain_mod_na)

par(mar = c(2, 2, 1, 1))
plot(domain_pred_na, col = clrs, breaks = brks, main = "Domain")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


# PRESENCE/BACKGROUND MODELLING: MAXENT ####

## MAXENT

### with altitude
?maxent  # 'maxent' function in 'dismo' package

maxent_mod_po <- maxent(x = layers_mod_presonly[[vars_sel]], p = presence_centroids)  # ERROR: it won't compute with the presence-only pixels, it needs to generate background points at non-presence pixels too, so it's not a presence-only method!
maxent_mod <- maxent(x = layers_mod[[vars_sel]], p = presence_centroids)
# you can also use the data frame rather than the raster maps + presence points:
#maxent_mod <- maxent(x = dat[ , vars_sel], p = dat$presence)

# if you get an error like "MaxEnt is missing or incompatible with your version of Java" (or something similar), try copying the 'maxent.jar' file provided in the 'practicals' folder to the folder that you get with:
#system.file("java", package = "dismo")
# then try the 'maxent' command above again. If it still doesn't work, go to Session - Restart R" and run the script up to this point again. And if it still doesn't work, don't worry - we'll also build Maxent models without Java below.

# let's look at the contents of the resulting object:
str(maxent_mod)
nrow(maxent_mod@presence)
nrow(maxent_mod@absence)  # it creates absences from the pixels without presence points! so it does use absence data

# compute and map maxent predictions:
maxent_pred <- predict(layers_mod, maxent_mod)
plot(maxent_pred, col = clrs, breaks = brks, main = "Maxent")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


### without altitude
maxent_mod_na <- maxent(x = layers_mod[[vars_sel2]], p = presence_centroids)

str(maxent_mod_na)

maxent_pred_na <- predict(layers_mod, maxent_mod_na)
plot(maxent_pred_na, col = clrs, breaks = brks, main = "Maxent - no altitude")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


## MAXNET

### with altitude
?maxnet  # 'maxnet' function in 'maxnet' package, which uses 'glmnet' rather than the Java version of Maxent

# let's first create a version of 'dat' with only the presence data:
dat_presonly <- dat[dat$presence == 1, ]

maxnet_mod_po <- maxnet(p = dat_presonly$presence, data = dat_presonly[ , vars_sel], maxnet.formula(p = dat_presonly$presence, data = dat_presonly[ , vars_sel], classes = "lq"))  # ERROR: also won't work with a presence-only dataframe; it needs non-presence rows too:

maxnet_mod <- maxnet(p = dat$presence, data = dat[ , vars_sel], maxnet.formula(p = dat$presence, data = dat[ , vars_sel]))

# compute and map maxnet predictions:
maxnet_pred <- predict(layers_mod, maxnet_mod, type = "cloglog")
plot(maxnet_pred, col = clrs, breaks = brks, main = "Maxnet")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


### without altitude
maxnet_mod_na <- maxnet(p = dat$presence, data = dat[ , vars_sel2], maxnet.formula(p = dat$presence, data = dat[ , vars_sel2]))

maxnet_pred_na <- predict(layers_mod, maxnet_mod_na, type = "cloglog")
plot(maxnet_pred_na, col = clrs, breaks = brks, main = "Maxnet - no altitude")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


# PRESENCE/ABSENCE MODELLING: GLM ####

## with altitude
?glm

# let's first make the formula required by the 'glm' function:
glm_form <- as.formula(paste("presence ~", paste(vars_sel, collapse = " + ")))
glm_form

glm_mod <- glm(formula = glm_form, family = binomial, data = dat)
summary(glm_mod, correlation = TRUE)

glm_pred <- predict(layers_mod, glm_mod, type = "response")  # "response" transforms the linear predictor with the appropriate link function, yielding results as probability values
plot(glm_pred, col = clrs, breaks = brks, main = "GLM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# note that the predictions of GLM and other presence-absence models are of presence PROBABILITY, which incorporates the effects of the predictor variables AND the baseline prevalence (proportion of presences) of the species in the model training data
# probabilities are generally low for relatively rare species, and generally higher for widespread species

# convert probability to prevalence-independent favourability:
glm_fav <- Fav(pred = glm_pred, sample.preval = prevalence(glm_mod$y))
plot(glm_fav, col = clrs, breaks = brks, main = "GLM favourability")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


## without altitude
glm_form2 <- as.formula(paste("presence ~", paste(vars_sel2, collapse = " + ")))
glm_form2

glm_mod_na <- glm(formula = glm_form2, family = binomial, data = dat)
summary(glm_mod_na, correlation = TRUE)

glm_pred_na <- predict(layers_mod, glm_mod_na, type = "response") 
plot(glm_pred_na, col = clrs, breaks = brks, main = "GLM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

glm_fav_na <- Fav(pred = glm_pred_na, sample.preval = prevalence(glm_mod_na$y))
plot(glm_fav_na, col = clrs, breaks = brks, main = "GLM favourability - no altitude")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


# PRESENCE/ABSENCE MODELLING: GAM ####

## with altitude
?gam

gam_form <- as.formula(paste("presence ~", paste0("s(", vars_sel, ")", collapse = "+")))  # GAM with smoothing splines ('s')
gam_form

gam_mod <- gam(formula = gam_form, family = binomial, data = dat)
summary(gam_mod)


gam_pred <- predict(layers_mod, gam_mod, type = "response")
plot(gam_pred, col = clrs, breaks = brks, main = "GAM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability:
gam_fav <- Fav(pred = gam_pred, sample.preval = prevalence(gam_mod$y))
plot(gam_fav, col = clrs, breaks = brks, main = "GAM favourability")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


## without altitude
gam_form2 <- as.formula(paste("presence ~", paste0("s(", vars_sel2, ")", collapse = "+")))  
gam_form2

gam_mod_na <- gam(formula = gam_form2, family = binomial, data = dat)
summary(gam_mod_na)


gam_pred_na <- predict(layers_mod, gam_mod_na, type = "response")
plot(gam_pred_na, col = clrs, breaks = brks, main = "GAM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

gam_fav_na <- Fav(pred = gam_pred_na, sample.preval = prevalence(gam_mod_na$y))
plot(gam_fav_na, col = clrs, breaks = brks, main = "GAM favourability - no altitude")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)



# PRESENCE/ABSENCE MODELLING: RANDOM FOREST (RF) ####

## with altitude
?randomForest
presence_fact <- as.factor(dat[ , "presence"])  # we need to convert the "presence" column from integer to factor (categorical variable), otherwise 'randomForest' would do regression as if for a continuous response variable (rather than binary presence/absence)

rf_form <- as.formula(paste("presence_fact ~", paste(vars_sel, collapse = "+")))
rf_form

rf_mod <- randomForest(formula = rf_form, data = dat, na.action = na.exclude)  # I've used the defaults, but you should read the help to explore important parameters if you're interested in using RF
rf_mod

rf_pred <- 1 - predict(layers_mod, rf_mod, type = "prob")
plot(rf_pred, col = clrs, breaks = brks, main = "RF")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)  # note that random forests tend to overfit to the observed records

# convert probability to prevalence-independent favourability:
rf_fav <- Fav(pred = rf_pred, sample.preval = prevalence(rf_mod$y))
plot(rf_fav, col = clrs, breaks = brks, main = "RF favourability")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# notice how RF overfits to stripes of presences missing from the main data source:
plot(rf_fav, xlim = c(-7, -4), ylim = c(36, 39), col = clrs, breaks = brks, main = "RF favourability")
plot(presence_centroids, pch = 20, col = "red", add = TRUE)


## without altitude
rf_form2 <- as.formula(paste("presence_fact ~", paste(vars_sel2, collapse = "+")))
rf_form2

rf_mod_na <- randomForest(formula = rf_form2, data = dat, na.action = na.exclude)  
rf_mod_na

rf_pred_na <- 1 - predict(layers_mod, rf_mod_na, type = "prob")
plot(rf_pred_na, col = clrs, breaks = brks, main = "RF")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)  

rf_fav_na <- Fav(pred = rf_pred_na, sample.preval = prevalence(rf_mod_na$y))
plot(rf_fav_na, col = clrs, breaks = brks, main = "RF favourability - no altitude")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# notice how RF overfits to stripes of presences missing from the main data source:
plot(rf_fav_na, xlim = c(-7, -4), ylim = c(36, 39), col = clrs, breaks = brks, main = "RF favourability - no altitude")
plot(presence_centroids, pch = 20, col = "red", add = TRUE)


# PRESENCE/ABSENCE MODELLING: BOOSTED REGRESSION TREES (BRT) ####

## with altitude
?gbm
gbm_form <- glm_form

gbm_mod <- gbm(formula = gbm_form, data = dat)  # I've used the defaults, but you should read the help to explore important parameters like 'shrinkage' (learning rate) and 'interaction.depth' (tree complexity) if you're interested in using GBM/BRT
gbm_mod

gbm_pred <- predict(layers_mod, gbm_mod, type = "response")
plot(gbm_pred, col = clrs, breaks = brks, main = "GBM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability:
gbm_fav <- Fav(pred = gbm_pred, sample.preval = prevalence(gbm_mod$data$y))
plot(gbm_fav, col = clrs, breaks = brks, main = "GBM favourability")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

## without altitude
gbm_form2 <- glm_form2

gbm_mod_na <- gbm(formula = gbm_form2, data = dat) 
gbm_mod_na

gbm_pred_na <- predict(layers_mod, gbm_mod_na, type = "response")
plot(gbm_pred_na, col = clrs, breaks = brks, main = "GBM")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)

# convert probability to prevalence-independent favourability:
gbm_fav_na <- Fav(pred = gbm_pred_na, sample.preval = prevalence(gbm_mod_na$data$y))
plot(gbm_fav_na, col = clrs, breaks = brks, main = "GBM favourability")
plot(presence_centroids, pch = ".", col = "red", add = TRUE)


# plot the raster map predictions together:

par(mar = c(2, 2, 2, 1), mfrow = c(4, 2))
plot(bioclim_pred, col = clrs, main = "Bioclim")
plot(domain_pred, col = clrs, main = "Domain")
plot(maxent_pred, col = clrs, main = "Maxent")
plot(maxnet_pred, col = clrs, main = "Maxnet")
plot(glm_fav, col = clrs, main = "GLM fav")
plot(gam_fav, col = clrs, main = "GAM fav")
plot(rf_fav, col = clrs, main = "RF fav")
plot(gbm_fav, col = clrs, main = "GBM fav")

par(mar = c(2, 2, 2, 1), mfrow = c(4, 2))
plot(bioclim_pred_na, col = clrs, main = "Bioclim - no altitude")
plot(domain_pred_na, col = clrs, main = "Domain - no altitude")
plot(maxent_pred_na, col = clrs, main = "Maxent - no altitude")
plot(maxnet_pred_na, col = clrs, main = "Maxnet - no altitude")
plot(glm_fav_na, col = clrs, main = "GLM fav - no altitude")
plot(gam_fav_na, col = clrs, main = "GAM fav - no altitude")
plot(rf_fav_na, col = clrs, main = "RF fav - no altitude")
plot(gbm_fav_na, col = clrs, main = "GBM fav - no altitude")



# APPLY PREDICTIONS TO THE DATA TABLE ####
# (variables must have exact same names as in the models)

dat$bioclim_pred <- predict(bioclim_mod, dat)
dat$domain_pred <- predict(domain_mod, dat)
dat$maxent_pred <- predict(maxent_mod, dat)
dat$maxnet_pred <- as.vector(predict(maxnet_mod, dat, type = "cloglog"))
dat$glm_pred <- predict(glm_mod, newdata = dat, type = "response")
dat$glm_fav <- Fav(pred = dat$glm_pred, sample.preval = prevalence(glm_mod$y))
dat$gam_pred <- predict(gam_mod, newdata = dat, type = "response")
dat$gam_fav <- Fav(pred = dat$gam_pred, sample.preval = prevalence(gam_mod$y))
dat$rf_pred <- predict(rf_mod, newdata = dat, type = "prob")[ , "1"]  # RF prediction includes probabilities for each observed class; we want the probability of the presence class ("1")
dat$rf_fav <- Fav(pred = dat$rf_pred, sample.preval = prevalence(rf_mod$y))
dat$gbm_pred <- predict(gbm_mod, newdata = dat, type = "response")
dat$gbm_fav <- Fav(pred = dat$gbm_pred, sample.preval = prevalence(gbm_mod$data$y))


dat$bioclim_pred_na <- predict(bioclim_mod_na, dat)
dat$domain_pred_na <- predict(domain_mod_na, dat)
dat$maxent_pred_na <- predict(maxent_mod_na, dat)
dat$maxnet_pred_na <- as.vector(predict(maxnet_mod_na, dat, type = "cloglog"))
dat$glm_pred_na <- predict(glm_mod_na, newdata = dat, type = "response")
dat$glm_fav_na <- Fav(pred = dat$glm_pred_na, sample.preval = prevalence(glm_mod_na$y))
dat$gam_pred_na <- predict(gam_mod_na, newdata = dat, type = "response")
dat$gam_fav_na <- Fav(pred = dat$gam_pred_na, sample.preval = prevalence(gam_mod_na$y))
dat$rf_pred_na <- predict(rf_mod_na, newdata = dat, type = "prob")[ , "1"]  
dat$rf_fav_na <- Fav(pred = dat$rf_pred_na, sample.preval = prevalence(rf_mod_na$y))
dat$gbm_pred_na <- predict(gbm_mod_na, newdata = dat, type = "response")
dat$gbm_fav_na <- Fav(pred = dat$gbm_pred_na, sample.preval = prevalence(gbm_mod_na$data$y))

head(dat)

pred_cols <- c("bioclim_pred", "domain_pred", "maxent_pred", "maxnet_pred", "glm_fav", "gam_fav", "rf_fav", "gbm_fav")
pred_cols_na <- c("bioclim_pred_na", "domain_pred_na", "maxent_pred_na", "maxnet_pred_na", "glm_fav_na", "gam_fav_na", "rf_fav_na", "gbm_fav_na")



# map predictions from the data table:

# first we'll need to convert 'dat' to a spatial object:
dat_spatial <- dat
names(dat_spatial)
coordinates(dat_spatial) <- dat_spatial[ , c("x", "y")]

# map all prediction columns in one window:
spplot(dat_spatial, zcol = pred_cols, col.regions = clrs, cex = 0.1)
spplot(dat_spatial, zcol = pred_cols_na, col.regions = clrs, cex = 0.1)

# map each prediction in turn:
spplot(dat_spatial, zcol = "bioclim_pred", cuts = brks, col.regions = clrs, cex = 0.5, main = "Bioclim")
spplot(dat_spatial, zcol = "domain_pred", cuts = brks, col.regions = clrs, cex = 0.5, main = "Domain")
spplot(dat_spatial, zcol = "maxent_pred", cuts = brks, col.regions = clrs, cex = 0.5, main = "Maxent")
spplot(dat_spatial, zcol = "maxnet_pred", cuts = brks, col.regions = clrs, cex = 0.5, main = "Maxnet")
spplot(dat_spatial, zcol = "glm_fav", cuts = brks, col.regions = clrs, cex = 0.5, main = "GLM fav")
spplot(dat_spatial, zcol = "gam_fav", cuts = brks, col.regions = clrs, cex = 0.5, main = "GAM fav")
spplot(dat_spatial, zcol = "rf_fav", cuts = brks, col.regions = clrs, cex = 0.5, main = "RF fav")
spplot(dat_spatial, zcol = "gbm_fav", cuts = brks, col.regions = clrs, cex = 0.5, main = "GBM fav")

spplot(dat_spatial, zcol = "bioclim_pred_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "Bioclim")
spplot(dat_spatial, zcol = "domain_pred_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "Domain")
spplot(dat_spatial, zcol = "maxent_pred_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "Maxent")
spplot(dat_spatial, zcol = "maxnet_pred_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "Maxnet")
spplot(dat_spatial, zcol = "glm_fav_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "GLM fav")
spplot(dat_spatial, zcol = "gam_fav_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "GAM fav")
spplot(dat_spatial, zcol = "rf_fav_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "RF fav")
spplot(dat_spatial, zcol = "gbm_fav_na", cuts = brks, col.regions = clrs, cex = 0.5, main = "GBM fav")


# COMPARE PREDICTIONS OF DIFFERENT MODELS ####

# pair-wise scatterplots of predictions:
pairs(dat[ , pred_cols], pch = 20, cex = 0.1)  # takes time!
pairs(dat[ , pred_cols_na], pch = 20, cex = 0.1)  # takes time!

# correlations among predictions:
pred_corrs <- cor(dat[ , pred_cols])
pred_corrs
min(pred_corrs, na.rm = TRUE)

pred_corrs_na <- cor(dat[ , pred_cols_na])
pred_corrs_na
min(pred_corrs_na, na.rm = TRUE)

# visual correlation matrix:
par(mfrow = c(1, 1))
corrplot(pred_corrs, method = "ellipse", type = "upper", addCoef.col = "wheat3", addCoefasPercent = TRUE)

par(mfrow = c(1, 1))
corrplot(pred_corrs_na, method = "ellipse", type = "upper", addCoef.col = "wheat3", addCoefasPercent = TRUE)


# save all current objects to disk (YOU WILL NEED THEM LATER!):
save.image(paste0("./data/niche_model/4modelling/", str_remove(myspecies, "Linum "), ".RData"))
