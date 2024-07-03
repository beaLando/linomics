library("phytools")


### Import tree files
bestML <- read.tree("./output/Angiosperm353/Astral/astral353.tre") # Best ML tree of the divaricate phylogeny
#BSreps <- read.tree("July20_cpDNA_BSreps.tre") # The corresponding BS replicates



##### (1) Determine the smallest BL of the best ML tree and BS replicates

### Best ML tree
min(bestML$edge.length, na.rm = TRUE)
# Result: 0
# Add 0.0000001 (br lengths have 6 decimals, so I add one more decimal to add very small number. I have to do it, because 0 gives always zero if multiplied)
# => Multiplying factor: x100 is enough to bring all values > 3e-6 for the best ML tree #3.5001383e-06 is the standardised value used by treePL when it finds short branches

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
bestMLx100$edge.length <- bestMLx100$edge.length+0.0000001
bestMLx100$edge.length <- bestMLx100$edge.length*100
write.tree(bestMLx100, "./output/Angiosperm353/Astral/astral353_modifiedBrL.tre")


### Multiply BL of the BS replicates
BSrepsx100 <- BSreps
for (i in c(1:1008)){ # Note that you usually obtain a few more BS replicates than you asked for from RAxML
  BSrepsx100[[i]]$edge.length <- BSrepsx100[[i]]$edge.length*100
}
write.tree(BSrepsx100, "July20_cpDNA_BSrepsx100.tre")