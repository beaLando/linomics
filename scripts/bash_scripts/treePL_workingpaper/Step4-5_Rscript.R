library("phytools", lib.loc="~/R/R-3.6.1/library")
setwd("path/to/your/working/directory")



### Import tree files
bestML <- read.tree("July20_cpDNA_noBS.testB") # Best ML tree of the divaricate phylogeny
BSreps <- read.tree("July20_cpDNA_BSreps.tre") # The corresponding BS replicates



##### (1) Determine the smallest BL of the best ML tree and BS replicates

### Best ML tree
min(bestML$edge.length)
# Result: 0.00000100.
# => Multiplying factor: x100 is enough to bring all values > 3e-5 for the best ML tree

### BS replicates
# You could create a function or a for loop that searches min(BSreps$edge.length) for each BS
# replicate then finds the lowest of all these minima, but I found it more efficient to open
# July20_cpDNA_BSreps.tre with a text editor and search for "0.000000" (one power of 10 under
# the minimum BL of the best ML tree) to check whether there was at least one BS replicates that
# had at least one BL value that would have needed a multiplying factor of 1000 or above instead
# of 100 to become longer than the minimum BL.



##### (2) Multiply the rates


### Multiply BL of the best ML tree
bestMLx100 <- bestML
bestMLx100$edge.length <- bestMLx100$edge.length*100
write.tree(bestMLx100, "July20_cpDNA_noBSx100.tre")


### Multiply BL of the BS replicates
BSrepsx100 <- BSreps
for (i in c(1:1008)){ # Note that you usually obtain a few more BS replicates than you asked for from RAxML
  BSrepsx100[[i]]$edge.length <- BSrepsx100[[i]]$edge.length*100
}
write.tree(BSrepsx100, "July20_cpDNA_BSrepsx100.tre")