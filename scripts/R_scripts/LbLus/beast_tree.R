# LIBRARIES
library(strap)
library(phytools)
library(phyloch)

# Compare beast output with IQtree
## BEAST
tree <- read.nexus("./output/plastome/beast/all_plastomes_aligned_trimmed_supermx_beast_17NOempty.(time).MCMC.trees")
plot(tree, cex = 0.5)

## IQtree
treP <- ape::read.tree(file = "./output/plastome/iqtree/without_84/formatted_consensus.nwk.txt")
plot(treP, cex = 0.5)

left <- ladderize(treP, right = T)
plot(left)
left2 <- drop.tip(left, "L88")
plot(left2)

right <- ladderize(tree, right = T)
plot(right)
right2 <- drop.tip(right, c("L88"))
plot(right2)

both.tre <- cophylo(left2, right2, rotate = TRUE)


dev.off()

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off() #only small differences in topology


# Plot BEAST output
tree$root.time <- max(nodeHeights(tree))
geoscalePhylo(ladderize(tree, right = F), 
              x.lim = c(0, 7), 
              cex.tip = 0.7, 
              cex.age = 1.3, 
              cex.ts = 1)

annot_tree <- phyloch::read.beast("./output/plastome/beast/all_plastomes_aligned_trimmed_supermx_beast_17NOempty.(time).MCMC.trees")

if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
  annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
} else {
  annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
}


annot_tree$root.time <- max(nodeHeights(annot_tree))
annot_tree$nodlabs <- as.character(1:94)

dt <- data.frame(nodlabs = annot_tree$nodlabs, 
                 min.ag = as.numeric(annot_tree$min_ages), 
                 max.ag = as.numeric(annot_tree$max_ages),
                 min.hpd = as.numeric(annot_tree$`height_95%_HPD_MIN`),
                 max.hpd = as.numeric(annot_tree$`height_95%_HPD_MAX`),
                 median.ag = as.numeric(annot_tree$`height_median`))

geoscalePhylo(ladderize(annot_tree, right = F),
              x.lim = c(0, 20), 
              cex.tip = 0.5, 
              cex.age = 1.3, 
              cex.ts = 0.7)

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(annot_tree) + 1):(annot_tree$Nnode + Ntip(annot_tree))) {
  lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(annot_tree)],
              T1$root.time - annot_tree$max_ages[i - Ntip(annot_tree)]),
        y = rep(T1$yy[i], 2), lwd = 1.5, lend = 0,
        col = make.transparent("blue", 0.4))
  }
