# LIBRARIES
library(ape)
library(phytools)
library(BAMMtools)
library(coda)

treepl_tree_matrix_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_g/supermx/supermx353_linum_sm0_primed_dated_ultraPAR.tre")

edata_matrix <- getEventData(treepl_tree_matrix_gtrg, 
                             "data/Angiosperm353/iqtrees/gtr_g/bamm_res/supermx_eventData.txt", 
                             burnin=0.1)


# CHECK MCMC CONVERGENCE
mcmcout_matrix <- read.csv("data/Angiosperm353/iqtrees/gtr_g/bamm_res/supermx_out.txt", header=T)
plot(mcmcout_matrix$logLik ~ mcmcout_matrix$generation)

# BURN-IN
burnstart_matrix <- floor(0.1 * nrow(mcmcout_matrix))
postburn_matrix <- mcmcout_matrix[burnstart_matrix:nrow(mcmcout_matrix), ]

effectiveSize(postburn_matrix$N_shifts) #ok >200
effectiveSize(postburn_matrix$logLik) #ok >200

# PROBABILITY of NUMBER of RATE SHIFTS
post_probs_matrix <- table(postburn_matrix$N_shifts) / nrow(postburn_matrix) #there is a higher probability than 0 rate shifts for 1 and 2 rate shifts
names(post_probs_matrix) 

shift_probs_matrix <- summary(edata_matrix) #another way to look at it

# PRIOR DISTRIBUTION
## if the bayesian factor for N shifts is above 12, that N is considered to improve the model relative to a 0 shift model
bfmat_matrix <- computeBayesFactors(mcmcout_matrix, expectedNumberOfShifts=1, burnin=0.1) #check first column
plotPrior(mcmcout_matrix, expectedNumberOfShifts=1) #in my case quite weak evidence

# PLOT MEAN PHYLORATE PLOT
plot.bammdata(edata_matrix, lwd=3, pal="temperature", legend = T)

# 95% MAXIMUM CREDIBLE SHIFTS and BEST CONFIGURATION
css_matrix <- credibleShiftSet(edata_matrix, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
summary(css_matrix)

css_matrix$number.distinct

plot.credibleshiftset(css_matrix)
dev.off()

best_matrix <- getBestShiftConfiguration(edata_matrix, expectedNumberOfShifts=1) #the number is given by previous step (in this case 1)
plot.bammdata(best_matrix, lwd = 2)
addBAMMshifts(best_matrix, cex=2.5)

# EXAMINE RANDOM SAMPLES FROM SHIFT CONFIGURATION
## You can see that, in each case, the same distinct shift configuration is being plotted, but a different sample. Hence, you should notice that the actual positions of shifts will move around.
dsc_matrix <- distinctShiftConfigurations(edata_matrix, expectedNumberOfShifts=1, threshold=5)
# Here is one random sample with the BEST shift configuration
plot.bammshifts(dsc_matrix, edata_matrix, rank=1, legend=F)
# Here is another (read the label text):
plot.bammshifts(dsc_matrix, edata_matrix, rank=1, legend=F)

# MARGINAL SHIFT PROBABILITIES
marg_probs_matrix <- marginalShiftProbsTree(edata_matrix)
plot.phylo(marg_probs_matrix)

# WHOLE CLADE AND SUB-CLADES EVOLUTIONARY RATES
## ALL
allrates_matrix <- getCladeRates(edata_matrix)
mean(allrates_matrix$lambda)
quantile(allrates_matrix$lambda, c(0.05, 0.95)) #0.12

plotRateThroughTime(edata_matrix, ratetype="speciation")

## LINUM
plot(treepl_tree_matrix_gtrg)
nodelabels()

LinumRates_matrix <- getCladeRates(edata_matrix, node=27)
mean(LinumRates_matrix$lambda)
quantile(LinumRates_matrix$lambda, c(0.05, 0.95)) #0.15

plotRateThroughTime(edata_matrix, node = 27, ratetype="speciation")

## BLUE
BLinumRates_matrix <- getCladeRates(edata_matrix, node=32)
mean(BLinumRates_matrix$lambda)
quantile(BLinumRates_matrix$lambda, c(0.05, 0.95)) #0.17

plotRateThroughTime(edata_matrix, node = 32, ratetype="speciation")

## YELLOW
YLinumRates_matrix <- getCladeRates(edata_matrix, node=28)
mean(YLinumRates_matrix$lambda)
quantile(YLinumRates_matrix$lambda, c(0.05, 0.95)) #0.16

plotRateThroughTime(edata_matrix, node = 28, ratetype="speciation")

## L. BIENNE NODE
LBLinumRates_matrix <- getCladeRates(edata_matrix, node=11)
mean(LBLinumRates_matrix$lambda)
quantile(LBLinumRates_matrix$lambda, c(0.05, 0.95)) #0.17

plotRateThroughTime(edata_matrix, node = 11, ratetype="speciation")

# COHORT ANALYSIS
cmat_matrix <- getCohortMatrix(edata_matrix)
cohorts(cmat_matrix, edata_matrix)
cst_matrix <- cumulativeShiftProbsTree(edata_matrix)
plot.phylo(cst_matrix)

edgecols_matrix  <- rep('black', length(treepl_tree_matrix_gtrg$edge.length))
is_highprobshift_matrix  <- cst_matrix $edge.length >= 0.95
edgecols_matrix[ is_highprobshift_matrix  ] <- "red"
plot.phylo(treepl_tree_matrix_gtrg, edge.color = edgecols_matrix )




















treepl_tree_astral_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_g/astral/astral353_linum_sm0_primed_dated_ultraPAR.tre")

edata_astral <- getEventData(treepl_tree_astral_gtrg, 
                             "data/Angiosperm353/iqtrees/gtr_g/bamm_res/astral_eventData.txt", 
                             burnin=0.1)


# CHECK MCMC CONVERGENCE
mcmcout_astral <- read.csv("data/Angiosperm353/iqtrees/gtr_g/bamm_res/astral_out.txt", header=T)
plot(mcmcout_astral$logLik ~ mcmcout_astral$generation)

# BURN-IN
burnstart_astral <- floor(0.1 * nrow(mcmcout_astral))
postburn_astral <- mcmcout_astral[burnstart_astral:nrow(mcmcout_astral), ]

effectiveSize(postburn_astral$N_shifts) #ok >200
effectiveSize(postburn_astral$logLik) #ok >200

# PROBABILITY of NUMBER of RATE SHIFTS
post_probs_astral <- table(postburn_astral$N_shifts) / nrow(postburn_astral) #there is a higher probability than 0 rate shifts for 1 and 2 rate shifts
names(post_probs_astral) 

shift_probs_astral <- summary(edata_astral) #another way to look at it

# PRIOR DISTRIBUTION
## if the bayesian factor for N shifts is above 12, that N is considered to improve the model relative to a 0 shift model
bfmat_astral <- computeBayesFactors(mcmcout_astral, expectedNumberOfShifts=1, burnin=0.1) #check first column
plotPrior(mcmcout_astral, expectedNumberOfShifts=1) #in my case quite weak evidence

# PLOT MEAN PHYLORATE PLOT
plot.bammdata(edata_astral, lwd=3, pal="temperature", legend = T)

# 95% MAXIMUM CREDIBLE SHIFTS and BEST CONFIGURATION
css_astral <- credibleShiftSet(edata_astral, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
summary(css_astral)

css_astral$number.distinct

plot.credibleshiftset(css_astral)
dev.off()

best_astral <- getBestShiftConfiguration(edata_astral, expectedNumberOfShifts=1) #the number is given by previous step (in this case 1)
plot.bammdata(best_astral, lwd = 2)
addBAMMshifts(best_astral, cex=2.5)

# EXAMINE RANDOM SAMPLES FROM SHIFT CONFIGURATION
## You can see that, in each case, the same distinct shift configuration is being plotted, but a different sample. Hence, you should notice that the actual positions of shifts will move around.
dsc_astral <- distinctShiftConfigurations(edata_astral, expectedNumberOfShifts=1, threshold=5)
# Here is one random sample with the BEST shift configuration
plot.bammshifts(dsc_astral, edata_astral, rank=1, legend=F)
# Here is another (read the label text):
plot.bammshifts(dsc_astral, edata_astral, rank=1, legend=F)

# MARGINAL SHIFT PROBABILITIES
marg_probs_astral <- marginalShiftProbsTree(edata_astral)
plot.phylo(marg_probs_astral)

# WHOLE CLADE AND SUB-CLADES EVOLUTIONARY RATES
## ALL
allrates_astral <- getCladeRates(edata_astral)
mean(allrates_astral$lambda)
quantile(allrates_astral$lambda, c(0.05, 0.95)) #0.12

plotRateThroughTime(edata_astral, ratetype="speciation")

## LINUM
plot(treepl_tree_astral_gtrg)
nodelabels()

LinumRates_astral <- getCladeRates(edata_astral, node=27)
mean(LinumRates_astral$lambda)
quantile(LinumRates_astral$lambda, c(0.05, 0.95)) #0.15

plotRateThroughTime(edata_astral, node = 27, ratetype="speciation")

## BLUE
BLinumRates_astral <- getCladeRates(edata_astral, node=32)
mean(BLinumRates_astral$lambda)
quantile(BLinumRates_astral$lambda, c(0.05, 0.95)) #0.17

plotRateThroughTime(edata_astral, node = 32, ratetype="speciation")

## YELLOW
YLinumRates_astral <- getCladeRates(edata_astral, node=28)
mean(YLinumRates_astral$lambda)
quantile(YLinumRates_astral$lambda, c(0.05, 0.95)) #0.16

plotRateThroughTime(edata_astral, node = 28, ratetype="speciation")

## L. BIENNE NODE
LBLinumRates_astral <- getCladeRates(edata_astral, node=38)
mean(LBLinumRates_astral$lambda)
quantile(LBLinumRates_astral$lambda, c(0.05, 0.95)) #0.17

plotRateThroughTime(edata_astral, node = 38, ratetype="speciation")



treepl_tree_astral_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_g/Astral/astral353_linum_sm0_primed_dated_ultraPAR.tre")

edata_astral <- getEventData(treepl_tree_astral_gtrg, 
                             "data/Angiosperm353/iqtrees/gtr_g/bamm_res/astral_eventData.txt", 
                             burnin=0.1)

plot.bammdata(edata_astral, lwd=3, method="polar", pal="temperature")


# COHORT ANALYSIS
cmat_astral <- getCohortMatrix(edata_astral)
cohorts(cmat_astral, edata_astral)
cst_astral <- cumulativeShiftProbsTree(edata_astral)
plot.phylo(cst_astral)

edgecols_astral  <- rep('black', length(treepl_tree_astral_gtrg$edge.length))
is_highprobshift_astral  <- cst_astral $edge.length >= 0.95
edgecols_astral[ is_highprobshift_astral  ] <- "red"
plot.phylo(treepl_tree_astral_gtrg, edge.color = edgecols_astral )