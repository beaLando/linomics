library(tidyverse)
library(reshape2)
library(ape)
library(poppr)
library(geiger)
library(phytools)
library(treeio)
library(ggplot2)
library(patchwork)
library(tidytext)
library(ggtree)

## nuclear tree
treN <- ape::read.tree(file = "./output/nuclear_lowDepth/genetic_distances_tree.nwk") #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
plot(treN, cex = 0.5)


## tfl1 tree #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
treF <- ape::read.tree(file = "./output/tfl1/tfl1_port_tree/consensus.rooted3.tre")

left <- ladderize(treN, right = T)
plot(left)


right <- ladderize(treF, right = T)
plot(right)

both.tre <- cophylo(left, right, rotate = TRUE)


dev.off()

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 2, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()



#Import population info
pop.info <- read.csv("./data/lblus_all_info.csv")

#Import MAFFT alignment trimmed in trimaAL
aln <- fasta2DNAbin("./data/Lb_TFL1/lb_port.mft.tal.fasta", quiet = FALSE, chunkSize = 10, snpOnly = FALSE)
aln <- DNAbin2genind(aln, pop = NULL, exp.char = c("a","t","g","c"), polyThres = 1/100)

ploidy(aln) <- 1

#aln.f <- missingno(aln, "loci")
aln.f <- missingno(aln, "geno", cutoff = 0.6)

names.fixed <- data.frame(inds = str_remove(indNames(aln.f), " 10003 bp")) %>%
  separate(inds, into = c("inds", "inds2"), sep = "L") %>% #View()
  mutate(
    inds = if_else(inds == "", "L", inds),
    inds2 = sprintf("%02d", as.numeric(inds2))
  ) %>%
  mutate(inds.fixed = paste0(inds, inds2)) %>%
  mutate(inds.fixed = str_remove(inds.fixed, "NA")) %>%
  as.data.frame()

indNames(aln.f) <- names.fixed$inds.fixed

as.data.frame(indNames(aln.f)) %>%
  filter(!`indNames(aln.f)` %in% pop.info$seq_code)

pop.info2 <- pop.info %>%
  filter(seq_code %in% indNames(aln.f)) %>% #str()
  mutate(
    region_purpose = if_else(spp == "usitatissimum", purpose, region)) %>%
  distinct(seq_code, .keep_all = TRUE) %>%
  arrange(factor(seq_code, levels = indNames(aln.f))) %>%
  as.data.frame() #%>% View()

aln.f@pop <- as.factor(pop.info2$region_purpose)
aln.f@strata <- pop.info2

# PRODUCE DAPC

#k = max
set.seed(435)
grp <- find.clusters(aln.f, max.n.clust = 20) #PCs = 30, K = 5

table(pop(aln.f), grp$grp)
table.value(table(pop(aln.f), grp$grp), col.lab = paste("k", 1:5))

set.seed(674)
dapc.aln <- dapc(aln.f, var.contrib = TRUE, scale = FALSE, pop =  grp$grp) #PCs = 30, D = 3

scatter(dapc.aln, xax = 1, yax = 2, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
scatter(dapc.aln, xax = 1, yax = 3, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
scatter(dapc.aln, xax = 2, yax = 3, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
assignplot(dapc.aln)

dapc.plot <- as.data.frame(dapc.aln$ind.coord)
dapc.plot$Group <- dapc.aln$grp
dapc.plot$pop <- aln.f$strata$pop
dapc.plot$ind <- aln.f$strata$seq_code
dapc.plot$spp <- aln.f$strata$spp
dapc.plot$lat <- aln.f$strata$lat
dapc.plot$purpose <- aln.f$strata$purpose
dapc.plot$reg <- aln.f$strata$region
dapc.plot$reg <- factor(dapc.plot$reg,levels(as.factor(dapc.plot$reg))[c(1, 9, 7, 6, 4, 5, 3, 8, 2)])
dapc.plot$hap <- aln.f$strata$LuTFL1
head(dapc.plot)

phylo.clust <- read.csv("./data/Lb_TFL1/phylogenetic_clusters_nuclear_plastome.csv")
str(phylo.clust)

dapc.plot <- dapc.plot %>%
  left_join(.,
            phylo.clust,
            by = c("ind" = "sample")) %>%
  as.data.frame()

#cols <- brewer.pal(n = nlevels(as.factor(dapc.res.m$dapc.group)), name = "Spectral")
#cols <- c("palevioletred4", "darkred", "brown1", "darkorange1", "darkgoldenrod2", "lightgoldenrod",   
          #"skyblue3", "steelblue4", "darkseagreen4", "gray31")
cols <- c("gray67", "white", "gray28", "lightgray", "black")

p1 <- ggplot(dapc.plot, aes(x = LD1, y = LD2, fill = Group, shape = group_nuclear))
p1 <- p1 + geom_point(size = 3, col = "black")
p1 <- p1 + theme_bw(base_size = 20) + theme(legend.position = "left")
p1 <- p1 + scale_shape_manual(values = c(22, 23, 24, 25)) #3, 11
p1 <- p1 + scale_colour_manual(values = cols)
p1 <- p1 + scale_fill_manual(values = cols)
p1

p2 <- ggplot(dapc.plot, aes(x = LD1, y = LD3, col = Group, fill = Group, shape = group_nuclear))
p2 <- p2 + geom_point(size = 3, col = "black")
p2 <- p2 + theme_bw(base_size = 15) + theme(legend.position = "left")
p2 <- p2 + scale_shape_manual(values = c(22, 23, 24, 25)) #3, 11
p2 <- p2 + scale_colour_manual(values = cols)
p2 <- p2 + scale_fill_manual(values = cols) + theme(legend.position = "none")
p2

p3 <- ggplot(dapc.plot, aes(x = LD2, y = LD3, col = Group, fill = Group, shape = group_nuclear))
p3 <- p3 + geom_point(size = 3, col = "black")
p3 <- p3 + theme_bw(base_size = 15) + theme(legend.position = "left")
p3 <- p3 + scale_shape_manual(values = c(22, 23, 24, 25)) #3, 11
p3 <- p3 + scale_colour_manual(values = cols)
p3 <- p3 + scale_fill_manual(values = cols) + theme(legend.position = "none")
p3

p1 + (p2/p3) + plot_layout(widths = c(2, 1))


ggplot(dapc.plot, aes(x = LD1, y = LD2, col = lat, fill = lat, shape = group_nuclear)) +
  geom_point(size = 3) +
  theme_bw(base_size = 20) + theme(legend.position = "left")

#dapc - group assignments
dapc.res <- as.data.frame(dapc.aln$posterior)
dapc.res$pop <- aln.f$strata$pop
dapc.res$ind <- aln.f$strata$seq_code
dapc.res$spp <- aln.f$strata$spp
dapc.res$reg <- aln.f$strata$region
dapc.res$lat <- as.character(aln.f$strata$lat)

dapc.res.m <-melt(dapc.res)
colnames(dapc.res.m)[6:7] <- c("dapc.group", "hpd")
dapc.res.m$reg <- factor(dapc.res.m$reg,levels(as.factor(dapc.res.m$reg))[c(1, 9, 7, 6, 4, 5, 3, 8, 2)])
dapc.res.m$p_i <- paste(dapc.res.m$pop, dapc.res.m$ind, sep="_")
dapc.res.m$lat <- as.numeric(dapc.res.m$lat)

dapc.res.m <- dapc.res.m %>%
  left_join(.,
            phylo.clust,
            by = c("ind" = "sample")) %>%
  mutate(group_nuclear = factor(group_nuclear,levels(as.factor(group_nuclear))[c(3, 4, 1, 2)])) %>%
  as.data.frame()

p4 <- ggplot(dapc.res.m, aes(x = reorder_within(p_i, lat, group_nuclear), y = hpd, fill = dapc.group)) + scale_x_reordered()
p4 <- p4 + geom_bar(stat='identity', alpha = 0.8, width = 0.8, col = "black")
p4 <- p4 + scale_fill_manual(values = cols) 
p4 <- p4 + facet_grid(~group_nuclear, scales = "free", space = "free")
p4 <- p4 + theme(panel.background = element_rect("white"), 
                 panel.grid.major = element_line("white"),
                 strip.background = element_rect(color="darkgrey", fill="white"),
                 strip.text.x = element_text(size = 10, color = "black"),
                 axis.text.x = element_text(size = 5, color = "black", angle = 45),
                 axis.title.x = element_blank(),
                 legend.position = "none")
p4