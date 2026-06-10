# R scripts for the course "Species distribution and ecological niche modelling in R"
# A. Marcia Barbosa (https://modtools.wordpress.com)


# DON'T FORGET:
# Session -> Set Working Directory -> To Source File Location


# LOAD THE NECESSARY PACKAGES IN THE CURRENT R SESSION ####
library(modEvA)
library(ecospat)
library(blockCV)
library(raster)
library(dismo)
library(gam)
library(maxnet)
library(randomForest)
library(gbm)
library(tidyverse)


# if you have (re)started a clean R session, load the models saved at the end of the previous practical:

spp <- c("linum bienne")

load(paste0("./data/niche_model/models_", "lblus", ".RData"))


# EVALUATE MODELS ON THE TRAINING DATA ####

names(dat.lb)

# area under the ROC Curve (AUC) ####

par(mfrow = c(2, 2), mar = c(3, 2, 2, 1))  # set plot margins and 3x2 plots in same window

with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxent_pred.lb), main = "Maxent"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxnet_pred.lb), main = "Maxnet"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(glm_pred.lb), main = "GLM"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(gam_pred.lb), main = "GAM"))

# area under the precision-recall Curve (AUC-PR) ####
par(mfrow = c(2, 2), mar = c(3, 2, 2, 1))  # set plot margins and 3x2 plots in same window

with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxent_pred.lb), curve = "PR", main = "Maxent"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxnet_pred.lb), curve = "PR", main = "Maxnet"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(glm_pred.lb), curve = "PR", main = "GLM"))
with(dat.lb, AUC(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(gam_pred.lb), curve = "PR", main = "GAM"))

# threshold-based classification metrics ####

?threshMeasures
classif_metrics <- c("CCR", "Sensitivity", "Specificity", "Precision", "Recall", "TSS", "kappa")

my_thresh <- "preval"

par(mfrow = c(2, 2), mar = c(5, 2, 2, 1))
with(dat.lb, threshMeasures(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxent_pred.lb), thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Maxent"))
with(dat.lb, threshMeasures(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxnet_pred.lb), thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Maxnet"))
with(dat.lb, threshMeasures(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(glm_pred.lb), thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "GLM"))
with(dat.lb, threshMeasures(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(gam_pred.lb), thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "GAM"))


# calculate and plot Miller calibration line:

par(mfrow = c(2, 2), mar = c(3, 2, 2, 1))
with(dat.lb, MillerCalib(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxent_pred.lb), main = "Maxent"))
with(dat.lb, MillerCalib(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(maxnet_pred.lb), main = "Maxnet"))
with(dat.lb, MillerCalib(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(glm_pred.lb), main = "GLM"))
with(dat.lb, MillerCalib(obs = presence_centroids.lb@data[,c(2:3)], pred = terra::rast(gam_pred.lb), main = "GAM"))


# COMPUTE BOYCE INDEX ####

# from raster variables + presence points:
# ecospat.boyce(obs = coordinates(presence_centroids.lb), fit = maxent_pred.lb)
# ecospat.boyce(obs = coordinates(presence_centroids.lb), fit = maxnet_pred.lb)
# ecospat.boyce(obs = coordinates(presence_centroids.lb), fit = glm_pred.lb)
# ecospat.boyce(obs = coordinates(presence_centroids.lb), fit = gam_pred.lb)

# from the data frame with observed and predicted values:
ecospat.boyce(obs = dat.lb[dat.lb$presence == 1, "maxent_pred"], fit = dat.lb[ , "maxent_pred"])
ecospat.boyce(obs = dat.lb[dat.lb$presence == 1, "maxnet_pred"], fit = dat.lb[ , "maxnet_pred"])
ecospat.boyce(obs = dat.lb[dat.lb$presence == 1, "glm_pred"], fit = dat.lb[ , "glm_pred"])
ecospat.boyce(obs = dat.lb[dat.lb$presence == 1, "gam_pred"], fit = dat.lb[ , "gam_pred"])

# BUT THIS JUST EVALUATES HOW THE MODELS FIT THE SAME DATA ON WHICH THEY WERE TRAINED
# YOU CAN SET ASIDE SOME DATA TO LEAVE OUT OF MODEL TRAINING AND USE FOR TESTING OUT-OF-SAMPLE PREDICTIVE CAPACITY
# BLOCK CROSS-VALIDATION (below) IS CURRENTLY THE MOST APPROPRIATE METHOD

# DIVIDE STUDY AREA INTO SPATIAL BLOCKS ####

names(dat.lb)
# convert 'dat' to a spatial object with its coordinate reference system:
dat_spatial.lb <- dat.lb
coordinates(dat_spatial.lb) <- dat.lb[ , c("x", "y")]
crs(dat_spatial.lb) <- crs(layers_mod.lb)
par(mfrow = c(1, 1))  # reset to 1 plot per window
plot(dat_spatial.lb)

# calculate the range of spatial autocorrelation in the modelling variables in the study area:
?spatialAutoRange
sarange.lb <- cv_spatial_autocor(subset(layers_mod.lb, vars_sel.lb))  # you can also use an additional argument, speciesData = subset(dat_spatial, presence == 1), but mind that this is not always recommended; see ?spatialAutoRange for more info
# click on the 'Plots' pane for the graphic results!
sarange.lb$range #648km

# blocks based on the autocorrelation range may be too large for the models to be able to capture the species-environment relationship adequately
# get spatial blocks of 150 km instead (or it could be size appropriate to ecology of species):
set.seed(123)  # set a particular seed of random numbers, so the next command alway provides the same set of blocks:
blocks.lb <- cv_spatial(x = dat_spatial.lb, column = "presence", r = layers_mod.lb[[1]], size = 325000, k = 5)  # you can use sarange$range as 'theRange' if you don't think it's too large
# you can also use an optional additional argument, species = "presence"; see ?spatialBlock

blocks.lb$folds_list
blocks.lb$folds_ids

# add the fold ID to the data frame:
dat.lb$foldID <- blocks.lb$folds_ids
head(dat.lb)

dat_spatial.lb$foldID <- blocks.lb$folds_ids
spplot(dat_spatial.lb, "foldID")  # each fold has a different colour


# COMPUTE MODELS AND GET PREDICTIONS LEAVING OUT EACH FOLD IN TURN ####
folds.lb <- sort(unique(dat.lb$foldID))
folds.lb

names(dat.lb)

for (f in folds.lb) {
  cat("modelling outside fold", f, "...\n")  # inform of progress
  dat_train.lb <- subset(dat.lb, foldID != f)

  glm_mod_fold.lb <- glm(formula = glm_form.lb, family = binomial, data = dat_train.lb)
  dat.lb[ , paste0("glm_fold", f, "_pred")] <- predict(glm_mod_fold.lb, dat.lb, type = "response")

  gam_mod_fold.lb <- gam(formula = gam_form.lb, family = binomial, data = dat_train.lb)
  dat.lb[ , paste0("gam_fold", f, "_pred")] <- predict(gam_mod_fold.lb, dat.lb, type = "response")

  maxnet_mod_fold.lb <- maxnet(p = dat_train.lb[ , "presence"], data = dat_train.lb[ , vars_sel.lb], f = maxnet.formula(dat_train.lb[ , "presence"], dat_train.lb[ , vars_sel.lb], classes = "lq"))  # using linear ('l') and quadratic ('q') features
  dat.lb[ , paste0("maxnet_fold", f, "_pred")] <- predict(maxnet_mod_fold.lb, dat.lb, type = "cloglog")

  maxent_mod_fold.lb <- maxent(p = dat_train.lb[ , "presence"], x = dat_train.lb[ , vars_sel.lb])  # using linear ('l') and quadratic ('q') features
  dat.lb[ , paste0("maxent_fold", f, "_pred")] <- predict(maxent_mod_fold.lb, dat.lb, type = "cloglog")
  
  gc()  # cleanup memory before next loop iteration
}  # end for s for f

# see the new predictions added to the the data frame:
head(dat.lb)


# EVALUATE EACH MODEL ON ITS VALIDATION FOLD ####
fold_cols.lb <- grep("_fold", names(dat.lb))
names(dat.lb)[fold_cols.lb]
models <- c("glm", "gam", "maxnet", "maxent")  # "bioclim", "domain", "maxent"
measures <- c("AUC", "AUPR", "TSS", "MCS")

# create an empty table to receive the cross-validation results:
crossval.lb <- as.data.frame(matrix(nrow = length(folds.lb), ncol = length(measures) * length(models)))
colnames(crossval.lb) <- c(outer(models, measures, FUN = paste, sep = "_"))
crossval.lb  # for now it's only filled with NAs

par(mfrow = c(5, 4))
for (m in models)  for (f in folds.lb) {
  fold_name <- paste0("fold", f)
  fold_col <- names(dat.lb)[grep(paste0(m, "_fold", f), names(dat.lb))]
  fold_dat <- subset(dat.lb, foldID == f)
  crossval.lb[f, paste(m, "AUC", sep = "_")] <- AUC(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = paste(m, "AUC"))
  crossval.lb[f, paste(m, "AUPR", sep = "_")] <- AUC(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], curve = "PR", simplif = TRUE, plot = TRUE, main = paste(m, "AUPR"))
  crossval.lb[f, paste(m, "TSS", sep = "_")] <- threshMeasures(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], thresh = my_thresh, measures = "TSS", simplif = TRUE, standardize = FALSE, main = paste(m, "TSS"))
  crossval.lb[f, paste(m, "MCS", sep = "_")] <- MillerCalib(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], main = paste(m, "Miller line"))$slope
}
# press the back arrow on the top left of your plotting window to see the fold evaluations for the different models

crossval.lb
crossval.lb <- crossval.lb[ , order(names(crossval.lb))]

# boxplots of cross-validation metrics:
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1))
boxplot(crossval.lb, col = rep(1:length(measures), times = length(models)), las = 2)
abline(h = 1, col = "green", lty = 2)  # remember Miller calibration slope (MCS) should ideally be close to 1 (not bigger = better)

# # we can also plot just the mean cross-validation performance for each model
# # (but note that the mean is quite limited information):
# crossval_means.lb <- sapply(crossval.lb, mean, na.rm = TRUE)
# barplot(crossval_means.lb, ylim = c(0, max(crossval.lb, na.rm = TRUE)), col = rep(1:length(models), each = length(measures)), las = 2)
# abline(h = 1, col = "darkgrey", lty = 2)

# possible thresholds for acceptable metric values:
# AUC: > 0.7 (Swets 1988, doi: 10.1126/science.3287615)
# TSS: > 0.4 (proportional to AUC 0.7, as TSS goes from -1 to 1)
# Miller slope: between 0.5 and 1.5 (Baquero et al. in review)

abline(h = 0.7, col = "black", lty = 2)
abline(h = 0.4, col = "blue", lty = 2)
abline(h = 0.5, col = "green", lty = 2)
abline(h = 1.5, col = "green", lty = 2)


# save the cross-validation results to a file on disk:
write.csv(crossval.lb, paste0("./data/niche_model/crossval_", spp[1], ".csv"), row.names = FALSE)
