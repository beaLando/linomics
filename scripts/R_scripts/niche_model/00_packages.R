# INSTALL NECESSARY R PACKAGES (IF THEY ARE NOT ALREADY INSTALLED) ####

# you may need recent versions of R and RStudio

if(!("uuid" %in% rownames(installed.packages())))
  install.packages("uuid")
if(!("wk" %in% rownames(installed.packages())))
  install.packages("wk")
if(!("httpcode" %in% rownames(installed.packages())))
  install.packages("httpcode")
if(!("rgbif" %in% rownames(installed.packages())))
  install.packages("rgbif")  
if(!("maps" %in% rownames(installed.packages())))
  install.packages("maps")
if(!("rgdal" %in% rownames(installed.packages())))
  install.packages("rgdal")
if(!("scrubr" %in% rownames(installed.packages()))) install.packages("scrubr")#this was not available for my R version so I installed development version with remotes::install_github("ropensci/scrubr")
if(!("raster" %in% rownames(installed.packages()))) install.packages("raster")
if(!("sdmpredictors" %in% rownames(installed.packages()))) install.packages("sdmpredictors")

if(!("dismo" %in% rownames(installed.packages()))) install.packages("dismo")
if(!("lattice" %in% rownames(installed.packages()))) install.packages("lattice")
if(!("maxnet" %in% rownames(installed.packages()))) install.packages("maxnet")
if(!("fuzzySim" %in% rownames(installed.packages()))) install.packages("fuzzySim", repos = "http://R-Forge.R-project.org")
if(!("gam" %in% rownames(installed.packages()))) install.packages("gam")
if(!("randomForest" %in% rownames(installed.packages()))) install.packages("randomForest")
if(!("gbm" %in% rownames(installed.packages()))) install.packages("gbm")
if(!("corrplot" %in% rownames(installed.packages()))) install.packages("corrplot")


if(!("blockCV" %in% rownames(installed.packages()))) install.packages("blockCV")
if(!("modEvA" %in% rownames(installed.packages()))) install.packages("modEvA", repos = "http://R-Forge.R-project.org")
if(!("ecospat" %in% rownames(installed.packages()))) install.packages("ecospat")
