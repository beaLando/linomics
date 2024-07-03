#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("treeio")
#BiocManager::install("ggtree")

library(ape)
library(geiger)
library(phytools)
library(treeio)
library(ggplot2)
library(ggtree)

tre.astral0 <- read.astral(file = "./output/Angiosperm353/Astral/gtr_g/astral353_rooted.tre")
tre.astral <- read.tree(file = "./output/Angiosperm353/Astral/gtr_g/astral353_rooted.tre")
tre.iq <- read.tree(file = "./output/Angiosperm353/iqtree/gtr_g/partitions_rooted/supermx353.contree")
plot(tre.astral, cex = 0.5)
plot(tre.iq, cex = 0.5)


#Get nodes of interest for dating
# syntax = mrca(phylogeny)["taxon1", "taxon2"]

mrca(tre.astral)["Phyllanthus_sp.","Linum_flavum"]
mrca(tre.iq)["Phyllanthus_sp.","Linum_flavum"]

mrca(tre.astral)["Linum_flavum","Linum_grandiflorum"]
mrca(tre.iq)["Linum_flavum","Linum_grandiflorum"]

######################################################################################
#CHRONOPL explanation

#If age.max = NULL (the default), it is assumed that age.min gives 
#exactly known ages. Otherwise, age.max and age.min must be of the 
#same length and give the intervals for each node. Some node may be 
#known exactly while the others are known within some bounds: the 
#values will be identical in both arguments for the former 
#(e.g., age.min = c(10, 5), age.max = c(10, 6), node = c(15, 18) 
#means that the age of node 15 is 10 units of time, and the age of 
#node 18 is between 5 and 6).'

######################################################################################

#examine difference between lambda values (lambda is smoothing in treePL)

#lambda = 0 
calib1a <- chronopl(tre.iq, lambda=0, age.min = c(33.9), 
                    age.max = c(37.2), node = c(28))
plot(calib1a, cex = 0.6)

#lambda = 1 
calib1b <- chronopl(tre.iq, lambda=1, age.min = c(33.9), 
                    age.max = c(37.2), node = c(28))
plot(calib1b, cex = 0.6)

#lambda = 1, phyllantoids_linoids node added
calib2b <- chronopl(tre.iq, lambda=1, age.min = c(94.8, 33.9), 
                    age.max = c(108.6, 37.2), node = c(27, 28))
plot(calib2b, cex = 0.6)


######################################################################################
#CHRONOS
tre.iq.rooted <- root(tre.iq, outgroup = "Euphorbia_mesembryanthemifolia")
mrca(tre.iq.rooted)["Linum_flavum","Linum_grandiflorum"]

tre.iq.timed0a <- chronos(tre.iq.rooted, lambda = 0, model = "relaxed", quiet = FALSE,
        calibration = makeChronosCalib(tre.iq.rooted, node = c(28), age.min = 33.9, age.max = 37.2, interactive = FALSE))

tre.iq.timed0b <- chronopl(tre.iq.rooted, lambda = 0, node = c(28), age.min = 33.9, age.max = 37.2)

#tre.iq.timed1 <- chronos(tre.iq.rooted, lambda = 1, model = "relaxed", quiet = FALSE,
                         #calibration = makeChronosCalib(tre.iq.rooted, node = c(28), age.min = 33.9, age.max = 37.2, interactive = FALSE))

plot(tre.iq.timed0b)
#plot(tre.iq.timed1)

p <- ggtree(tre.iq.timed0b) + geom_nodelab() + geom_nodepoint() + geom_tippoint() + geom_tiplab() + theme_tree2() + scale_x_continuous(labels = abs)
revts(p)

write.tree(tre.iq.timed0b, file = "./output/Angiosperm353/treePL/R_alternative_treePL.newick")
