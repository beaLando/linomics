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

myspecies <- c("Linum bienne")

load(paste0("./data/niche_model/4modelling/", str_remove(myspecies, "Linum "), ".RData"))


# EVALUATE MODELS ON THE TRAINING DATA ####

names(dat)


# area under the ROC Curve (AUC) ####

par(mfrow = c(4, 2), mar = c(3, 2, 2, 1))  # set plot margins and 3x2 plots in same window
?AUC
with(dat, AUC(obs = presence, pred = bioclim_pred, main = "Bioclim"))
with(dat, AUC(obs = presence, pred = domain_pred, main = "Domain"))
with(dat, AUC(obs = presence, pred = maxent_pred, main = "Maxent"))
with(dat, AUC(obs = presence, pred = maxnet_pred, main = "Maxnet"))
with(dat, AUC(obs = presence, pred = glm_pred, main = "GLM"))
with(dat, AUC(obs = presence, pred = gam_pred, main = "GAM"))
with(dat, AUC(obs = presence, pred = rf_pred, main = "RF"))
with(dat, AUC(obs = presence, pred = gbm_pred, main = "GBM"))

with(dat, AUC(obs = presence, pred = bioclim_pred_na, main = "Bioclim"))
with(dat, AUC(obs = presence, pred = domain_pred_na, main = "Domain"))
with(dat, AUC(obs = presence, pred = maxent_pred_na, main = "Maxent"))
with(dat, AUC(obs = presence, pred = maxnet_pred_na, main = "Maxnet"))
with(dat, AUC(obs = presence, pred = glm_pred_na, main = "GLM"))
with(dat, AUC(obs = presence, pred = gam_pred_na, main = "GAM"))
with(dat, AUC(obs = presence, pred = rf_pred_na, main = "RF"))
with(dat, AUC(obs = presence, pred = gbm_pred_na, main = "GBM"))


# area under the precision-recall Curve (AUC-PR) ####
par(mfrow = c(4, 2), mar = c(3, 2, 2, 1))
with(dat, AUC(obs = presence, pred = bioclim_pred, curve = "PR", main = "Bioclim"))
with(dat, AUC(obs = presence, pred = domain_pred, curve = "PR", main = "Domain"))
with(dat, AUC(obs = presence, pred = maxent_pred, curve = "PR", main = "Maxent"))
with(dat, AUC(obs = presence, pred = maxnet_pred, curve = "PR", main = "Maxnet"))
with(dat, AUC(obs = presence, pred = glm_pred, curve = "PR", main = "GLM"))
with(dat, AUC(obs = presence, pred = gam_pred, curve = "PR", main = "GAM"))
with(dat, AUC(obs = presence, pred = rf_pred, curve = "PR", main = "RF"))
with(dat, AUC(obs = presence, pred = gbm_pred, curve = "PR", main = "GBM"))

with(dat, AUC(obs = presence, pred = bioclim_pred_na, curve = "PR", main = "Bioclim"))
with(dat, AUC(obs = presence, pred = domain_pred_na, curve = "PR", main = "Domain"))
with(dat, AUC(obs = presence, pred = maxent_pred_na, curve = "PR", main = "Maxent"))
with(dat, AUC(obs = presence, pred = maxnet_pred_na, curve = "PR", main = "Maxnet"))
with(dat, AUC(obs = presence, pred = glm_pred_na, curve = "PR", main = "GLM"))
with(dat, AUC(obs = presence, pred = gam_pred_na, curve = "PR", main = "GAM"))
with(dat, AUC(obs = presence, pred = rf_pred_na, curve = "PR", main = "RF"))
with(dat, AUC(obs = presence, pred = gbm_pred_na, curve = "PR", main = "GBM"))


# threshold-based classification metrics ####

?threshMeasures
classif_metrics <- c("CCR", "Sensitivity", "Specificity", "Precision", "Recall", "TSS", "kappa")

# you can choose a threshold that optimizes classification performance for a particular measure (e.g. TSS) for each method
?optiThresh
my_thresh <- with(dat, optiThresh(obs = presence, pred = gbm_pred_na, interval = 0.001, measures = "TSS", optimize = "each"))
my_thresh
my_thresh <- my_thresh$optimals.each$threshold
my_thresh

# or you can instead use species prevalence as the threshold for everyone:
my_thresh <- "preval"

par(mfrow = c(4, 2), mar = c(5, 2, 2, 1))
with(dat, threshMeasures(obs = presence, pred = bioclim_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Bioclim"))
with(dat, threshMeasures(obs = presence, pred = domain_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Domain"))
with(dat, threshMeasures(obs = presence, pred = maxent_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Maxent"))
with(dat, threshMeasures(obs = presence, pred = maxnet_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "Maxnet"))
with(dat, threshMeasures(obs = presence, pred = glm_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "GLM"))
with(dat, threshMeasures(obs = presence, pred = gam_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "GAM"))
with(dat, threshMeasures(obs = presence, pred = rf_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "RF"))
with(dat, threshMeasures(obs = presence, pred = gbm_pred, thresh = my_thresh, measures = classif_metrics, ylim = c(0, 1), main = "GBM"))

# calculate and plot Miller calibration line:
?MillerCalib
par(mfrow = c(4, 2), mar = c(3, 2, 2, 1))
with(dat, MillerCalib(obs = presence, pred = bioclim_pred, main = "Bioclim"))
with(dat, MillerCalib(obs = presence, pred = domain_pred, main = "Domain"))
with(dat, MillerCalib(obs = presence, pred = maxent_pred, main = "Maxent"))
with(dat, MillerCalib(obs = presence, pred = maxnet_pred, main = "Maxnet"))
with(dat, MillerCalib(obs = presence, pred = glm_pred, main = "GLM"))
with(dat, MillerCalib(obs = presence, pred = gam_pred, main = "GAM"))
with(dat, MillerCalib(obs = presence, pred = rf_pred, main = "RF"))
with(dat, MillerCalib(obs = presence, pred = gbm_pred, main = "GBM"))


# COMPUTE BOYCE INDEX ####

?ecospat.boyce

# from raster variables + presence points:
ecospat.boyce(obs = coordinates(presence_centroids), fit = bioclim_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = domain_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = maxent_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = maxnet_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = glm_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = gam_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = rf_pred)
ecospat.boyce(obs = coordinates(presence_centroids), fit = gbm_pred)

# from the data frame with observed and predicted values:
ecospat.boyce(obs = dat[dat$presence == 1, "bioclim_pred"], fit = dat[ , "bioclim_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "domain_pred"], fit = dat[ , "domain_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "maxent_pred"], fit = dat[ , "maxent_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "maxnet_pred"], fit = dat[ , "maxnet_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "glm_pred"], fit = dat[ , "glm_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "gam_pred"], fit = dat[ , "gam_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "rf_pred"], fit = dat[ , "rf_pred"])
ecospat.boyce(obs = dat[dat$presence == 1, "gbm_pred"], fit = dat[ , "gbm_pred"])


# BUT THIS JUST EVALUATES HOW THE MODELS FIT THE SAME DATA ON WHICH THEY WERE TRAINED
# YOU CAN SET ASIDE SOME DATA TO LEAVE OUT OF MODEL TRAINING AND USE FOR TESTING OUT-OF-SAMPLE PREDICTIVE CAPACITY
# BLOCK CROSS-VALIDATION (below) IS CURRENTLY THE MOST APPROPRIATE METHOD
# THE R CODE IS A BIT COMPLEX, BUT DON'T WORRY IF YOU DON'T UNDERSTAND ALL OF IT
# it's basically the same things we've done so far, but looped over several folds of cross-validation blocks


# DIVIDE STUDY AREA INTO SPATIAL BLOCKS ####

?spatialBlock

names(dat)
# convert 'dat' to a spatial object with its coordinate reference system:
dat_spatial <- dat
coordinates(dat_spatial) <- dat[ , c("x", "y")]
crs(dat_spatial) <- crs(layers_mod)
par(mfrow = c(1, 1))  # reset to 1 plot per window
plot(dat_spatial)

# calculate the range of spatial autocorrelation in the modelling variables in the study area:
?spatialAutoRange
sarange <- spatialAutoRange(subset(layers_mod, vars_sel2))  # you can also use an additional argument, speciesData = subset(dat_spatial, presence == 1), but mind that this is not always recommended; see ?spatialAutoRange for more info
# click on the 'Plots' pane for the graphic results!
sarange$range

# blocks based on the autocorrelation range may be too large for the models to be able to capture the species-environment relationship adequately
# get spatial blocks of 150 km instead (or it could be size appropriate to ecology of species):
set.seed(123)  # set a particular seed of random numbers, so the next command alway provides the same set of blocks:
blocks <- spatialBlock(speciesData = dat_spatial, rasterLayer = layers_mod[[1]], theRange = 150000, k = 5)  # you can use sarange$range as 'theRange' if you don't think it's too large
# you can also use an optional additional argument, species = "presence"; see ?spatialBlock

blocks$folds
blocks$foldID

# add the fold ID to the data frame:
dat$foldID <- blocks$foldID
head(dat)

dat_spatial$foldID <- blocks$foldID
spplot(dat_spatial, "foldID")  # each fold has a different colour


# COMPUTE MODELS AND GET PREDICTIONS LEAVING OUT EACH FOLD IN TURN ####

folds <- sort(unique(dat$foldID))
folds

names(dat)

for (f in folds) {
  cat("modelling outside fold", f, "...\n")  # inform of progress
  dat_train <- subset(dat, foldID != f)

  # dat_train_presonly <- subset(dat_train, presence == 1)
  #
  # bioclim_mod_fold <- bioclim(x = as.matrix(dat_train_presonly[ , vars_sel]), p = as.matrix(dat_train_presonly[ , c("x", "y")]))
  # dat[ , paste0("bioclim_fold", f, "_pred")] <- predict(bioclim_mod_fold, newdata = dat)
  #
  # domain_mod_fold <- domain(x = as.matrix(dat_train_presonly[ , vars_sel]), p = as.matrix(dat_train_presonly[ , c("x", "y")]))
  # dat[ , paste0("domain_fold", f, "_pred")] <- predict(domain_mod_fold, newdata = dat)

  glm_mod_fold <- glm(formula = glm_form, family = binomial, data = dat_train)
  dat[ , paste0("glm_fold", f, "_pred")] <- predict(glm_mod_fold, dat, type = "response")

  gam_mod_fold <- gam(formula = gam_form, family = binomial, data = dat_train)
  dat[ , paste0("gam_fold", f, "_pred")] <- predict(gam_mod_fold, dat, type = "response")

  maxent_mod_fold <- maxent(p = dat_train[ , "presence"], x = dat_train[ , vars_sel])
  dat[ , paste0("maxent_fold", f, "_pred")] <- predict(maxent_mod_fold, dat)
  
  maxnet_mod_fold <- maxnet(p = dat_train[ , "presence"], data = dat_train[ , vars_sel], f = maxnet.formula(dat_train[ , "presence"], dat_train[ , vars_sel], classes = "lq"))  # using linear ('l') and quadratic ('q') features
  dat[ , paste0("maxnet_fold", f, "_pred")] <- predict(maxnet_mod_fold, dat, type = "cloglog")

  presence_fact_fold <- as.factor(dat_train[ , "presence"])  # we need to convert presence/absence from integer to factor (categorical variable), otherwise 'randomForest' would do regression as if for a continuous response variable (rather than binary presence-absence)
  rf_form_fold <- as.formula(paste("presence_fact_fold ~", paste(vars_sel, collapse = "+")))
  rf_mod_fold <- randomForest(formula = rf_form_fold, data = dat_train, na.action = na.exclude)
  dat[ , paste0("rf_fold", f, "_pred")] <- predict(rf_mod_fold, newdata = dat, type = "prob")[ , "1"]  # RF prediction includes probabilities for each observed class; we want the probability of the presence class ("1")

  gbm_mod_fold <- gbm(gbm_form, distribution = "bernoulli", data = dat_train)
  dat[ , paste0("gbm_fold", f, "_pred")] <- predict(gbm_mod_fold, newdata = dat, type = "response")

  gc()  # cleanup memory before next loop iteration
}  # end for s for f

# see the new predictions added to the the data frame:
head(dat)


# EVALUATE EACH MODEL ON ITS VALIDATION FOLD ####

fold_cols <- grep("_fold", names(dat))
names(dat)[fold_cols]
models <- c("glm", "gam", "maxent", "maxnet", "rf", "gbm")  # "bioclim", "domain"
measures <- c("AUC", "AUPR", "TSS", "MCS")

# create an empty table to receive the cross-validation results:
crossval <- as.data.frame(matrix(nrow = length(folds), ncol = length(measures) * length(models)))
colnames(crossval) <- c(outer(models, measures, FUN = paste, sep = "_"))
crossval  # for now it's only filled with NAs

par(mfrow = c(5, 3), mar = c(1.5, 2, 1, 1))
for (m in models)  for (f in folds) {
  fold_name <- paste0("fold", f)
  fold_col <- names(dat)[grep(paste0(m, "_fold", f), names(dat))]
  fold_dat <- subset(dat, foldID == f)
  crossval[f, paste(m, "AUC", sep = "_")] <- AUC(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = paste(m, "AUC"))
  crossval[f, paste(m, "AUPR", sep = "_")] <- AUC(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], curve = "PR", simplif = TRUE, plot = TRUE, main = paste(m, "AUPR"))
  crossval[f, paste(m, "TSS", sep = "_")] <- threshMeasures(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], thresh = my_thresh, measures = "TSS", simplif = TRUE, standardize = FALSE, main = paste(m, "TSS"))
  crossval[f, paste(m, "MCS", sep = "_")] <- MillerCalib(obs = fold_dat[ , "presence"], pred = fold_dat[ , fold_col], main = paste(m, "Miller line"))$slope
}
# press the back arrow on the top left of your plotting window to see the fold evaluations for the different models

crossval
crossval <- crossval[ , order(names(crossval))]

# boxplots of cross-validation metrics:
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1))
boxplot(crossval, col = rep(1:length(models), each = length(measures)), las = 2)
abline(h = 1, col = "darkgrey", lty = 2)  # remember Miller calibration slope (MCS) should ideally be close to 1 (not bigger = better)

# we can also plot just the mean cross-validation performance for each model
# (but note that the mean is quite limited information):
crossval_means <- sapply(crossval, mean, na.rm = TRUE)
barplot(crossval_means, ylim = c(0, max(crossval, na.rm = TRUE)), col = rep(1:length(models), each = length(measures)), las = 2)
abline(h = 1, col = "darkgrey", lty = 2)

# possible thresholds for acceptable metric values:
# AUC: > 0.7 (Swets 1988, doi: 10.1126/science.3287615)
# TSS: > 0.4 (proportional to AUC 0.7, as TSS goes from -1 to 1)
# Miller slope: between 0.5 and 1.5 (Baquero et al. in review)

abline(h = 0.7, col = "darkred", lty = 2)
abline(h = 0.4, col = "darkblue", lty = 2)
abline(h = 0.5, col = "orange", lty = 2)
abline(h = 1.5, col = "orange", lty = 2)


# save the cross-validation results to a file on disk:
write.csv(crossval, paste0("./output/niche_model/crossval_", str_remove(myspecies, "Linum "), ".csv"), row.names = FALSE)
