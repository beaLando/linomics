# LIBRARIES
library(ape)
library(phytools)
library(BAMMtools)

# TREES
treepl_tree_matrix_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_g/supermx/supermx353_linum_sm0_primed_dated.tre")
treepl_tree_astral_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_g/Astral/astral353_linum_sm0_primed_dated.tre")

plot(treepl_tree_matrix_gtrg)

# CHECKs
## MATRIX
#treepl_tree_matrix_gtrg <- drop.tip(treepl_tree_matrix_gtrg, c("euphorbia"))

plot(treepl_tree_matrix_gtrg)
treepl_tree_matrix_gtrg <- phytools::force.ultrametric(treepl_tree_matrix_gtrg) #i correct it with phytools
is.ultrametric(treepl_tree_matrix_gtrg) #ok
is.binary(treepl_tree_matrix_gtrg) #ok
min(treepl_tree_matrix_gtrg$edge.length) #>0 #ok

## ASTRAL
#treepl_tree_astral_gtrg <- drop.tip(treepl_tree_astral_gtrg, c("euphorbia"))
plot(treepl_tree_astral_gtrg)
treepl_tree_astral_gtrg <- phytools::force.ultrametric(treepl_tree_astral_gtrg) #i correct it with phytools
is.ultrametric(treepl_tree_astral_gtrg) #now ok
is.binary(treepl_tree_astral_gtrg) #ok
min(treepl_tree_astral_gtrg$edge.length) #>0 but very small

write.tree(treepl_tree_astral_gtrg, "./data/Angiosperm353/iqtrees/gtr_g/supermx/supermx353_linum_sm0_primed_dated_ultraPAR.tre")
write.tree(treepl_tree_astral_gtrg, "./data/Angiosperm353/iqtrees/gtr_g/Astral/astral353_linum_sm0_primed_dated_ultraPAR.tre")

# SET PRIORS for BOTH TREES
setBAMMpriors(read.tree("./data/Angiosperm353/iqtrees/gtr_g/supermx/supermx353_linum_sm0_primed_dated_ultraPAR.tre"),
              outfile = "data/Angiosperm353/iqtrees/gtr_g/bamm_res/supermx_priors.txt")

setBAMMpriors(read.tree("./data/Angiosperm353/iqtrees/gtr_g/Astral/astral353_linum_sm0_primed_dated_ultraPAR.tre"),
              outfile = "data/Angiosperm353/iqtrees/gtr_g/bamm_res/astral_priors.txt")
