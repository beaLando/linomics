library(vcfR)
library(poppr)
library(reshape2)
library(ggplot2)
library(popgraph) #install.packages("devtools"), install_github("dyerlab/popgraph")
library(RColorBrewer)
library(igraph)
library(tidyverse)


#Import population info
pop.info <- read.csv("./data/lblus_all_info.csv")

#Import MAFFT alignment trimmed in trimaAL
aln <- fasta2DNAbin("./data/Lb_plastomes_raw/Consensus/all_plastomes_aligned_trimmed.fa", quiet = FALSE, chunkSize = 10, snpOnly = FALSE)
aln <- DNAbin2genind(aln, pop = NULL, exp.char = c("a","t","g","c"), polyThres = 1/100)

ploidy(aln) <- 1

aln.f <- missingno(aln, "loci")
aln.f <- missingno(aln.f, "geno")

aln.f <- aln.f[indNames(aln.f) != "L84 156703 bp"]
aln.f <- aln.f[indNames(aln.f) != "L88 156703 bp"]

names.fixed <- str_remove(indNames(aln.f), " 156703 bp")
indNames(aln.f) <- names.fixed

pop.info2 <- pop.info %>%
  filter(seq_code %in% indNames(aln.f)) %>%
  as.data.frame()

aln.f@pop <- as.factor(pop.info2$region)
aln.f@strata <- pop.info2

#dapc

#k = max
set.seed(435)
grp <- find.clusters(aln.f, max.n.clust = 40) #PCs = 15, K = 7

table(pop(aln.f), grp$grp)
table.value(table(pop(aln.f), grp$grp), col.lab = paste("k", 1:7))

set.seed(674)
dapc.aln <- dapc(aln.f, var.contrib = TRUE, scale = FALSE, pop =  grp$grp) #PCs = 15, D = 2

scatter(dapc.aln, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
assignplot(dapc.aln)

dapc.plot <- as.data.frame(dapc.aln$ind.coord)
dapc.plot$Group <- dapc.aln$grp
dapc.plot$pop <- aln.f$strata$pop
dapc.plot$ind <- aln.f$strata$seq_code
dapc.plot$spp <- aln.f$strata$spp
dapc.plot$reg <- aln.f$strata$region
dapc.plot$reg <- factor(dapc.plot$reg,levels(as.factor(dapc.plot$reg))[c(1, 10, 8, 7, 5, 6, 4, 9, 3, 2)])
head(dapc.plot)

#cols <- brewer.pal(n = nlevels(as.factor(dapc.res.m$dapc.group)), name = "Spectral")
#cols <- c("palevioletred4", "darkred", "brown1", "darkorange1", "darkgoldenrod2", "lightgoldenrod",   
            #"lightcyan", "slategray1", "skyblue3", "steelblue4", "darkslategray", "darkseagreen4", "gray31")
cols <- c("skyblue3", "gray", "darkred", "darkseagreen4", "brown1", "darkgoldenrod2", "steelblue4")

p1 <- ggplot(dapc.plot, aes(x = LD1, y = LD2, col = Group, fill = Group, shape = reg))
p1 <- p1 + geom_point(size = 4)
p1 <- p1 + theme_bw(base_size = 20) + theme(legend.position = "left")
p1 <- p1 + scale_shape_manual(values = c(8, 24, 4, 21, 1, 23, 22, 25, 3, 11))
p1 <- p1 + scale_colour_manual(values = cols)
p1 <- p1 + scale_fill_manual(values = cols)
p1

i1 <- p1 + coord_cartesian(xlim = c(-84, -77.5), ylim = c(3, 7.5)) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
i2 <- p1 + coord_cartesian(xlim = c(54, 67.5), ylim = c(-30, -3)) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
i3 <- p1 + coord_cartesian(xlim = c(112.5, 120), ylim = c(82, 92)) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

p1 + (i3/i2/i1)


#dapc - group assignments
dapc.res <- as.data.frame(dapc.aln$posterior)
dapc.res$pop <- aln.f$strata$pop
dapc.res$ind <- aln.f$strata$seq_code
dapc.res$spp <- aln.f$strata$spp
dapc.res$reg <- aln.f$strata$region

dapc.res.m <-melt(dapc.res)
colnames(dapc.res.m)[5:6] <- c("dapc.group", "hpd")
dapc.res.m$reg <- factor(dapc.res.m$reg,levels(as.factor(dapc.res.m$reg))[c(1, 10, 8, 7, 5, 6, 4, 9, 3, 2)])
dapc.res.m$p_i <- paste(dapc.res.m$pop, dapc.res.m$ind, sep="_")

p2 <- ggplot(dapc.res.m, aes(x = p_i, y = hpd, fill = dapc.group))
p2 <- p2 + geom_bar(stat='identity', alpha = 0.8)
p2 <- p2 + scale_fill_manual(values = cols) 
p2 <- p2 + facet_grid(~reg, scales = "free", space = "free")
p2 <- p2 + theme(panel.background = element_rect("white"), 
               panel.grid.major = element_line("darkgrey"),
               strip.background = element_rect(color="darkgrey", fill="white"),
               strip.text.x = element_text(size = 10, color = "black"),
               axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "none")
p2

((p1 + (i3/i2/i1) + plot_layout(widths = c(2, 1))) / p2) + plot_layout(heights = c(5, 1), guides = "collect")


#loadings DAPC
set.seed(4)
contrib1 <- loadingplot(dapc.aln$var.contr, axis = 1, thres = 0.015, lab.jitter = 1)

set.seed(6)
contrib2 <- loadingplot(dapc.aln$var.contr, axis = 2, thres = 0.015, lab.jitter = 1)



#dapc

#k = 4
set.seed(435)
grp4 <- find.clusters(aln.f, max.n.clust = 40) #PCs = 15, K = 4

table(pop(aln.f), grp$grp)
table.value(table(pop(aln.f), grp4$grp), col.lab = paste("k", 1:4))

set.seed(674)
dapc4.aln <- dapc(aln.f, var.contrib = TRUE, scale = FALSE, pop =  grp4$grp) #PCs = 15, D = 2

scatter(dapc4.aln, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
assignplot(dapc4.aln)

dapc4.plot <- as.data.frame(dapc4.aln$ind.coord)
dapc4.plot$Group <- dapc4.aln$grp
dapc4.plot$pop <- aln.f$strata$pop
dapc4.plot$ind <- aln.f$strata$seq_code
dapc4.plot$spp <- aln.f$strata$spp
dapc4.plot$reg <- aln.f$strata$region
dapc4.plot$reg <- factor(dapc4.plot$reg,levels(as.factor(dapc4.plot$reg))[c(1, 10, 8, 7, 5, 6, 4, 9, 3, 2)])
head(dapc4.plot)

#cols <- brewer.pal(n = nlevels(as.factor(dapc4.res.m$dapc4.group)), name = "Spectral")
#cols <- c("palevioletred4", "darkred", "brown1", "darkorange1", "darkgoldenrod2", "lightgoldenrod",   
#"lightcyan", "slategray1", "skyblue3", "steelblue4", "darkslategray", "darkseagreen4", "gray31")
cols <- c("darkgoldenrod2", "darkseagreen4", "brown1", "steelblue4")


p3 <- ggplot(dapc4.plot, aes(x = LD1, y = LD2, col = Group, fill = Group, shape = reg))
p3 <- p3 + geom_point(size = 4)
p3 <- p3 + theme_bw(base_size = 20) + theme(legend.position = "left")
p3 <- p3 + scale_shape_manual(values = c(8, 24, 4, 21, 1, 23, 22, 25, 3, 11))
p3 <- p3 + scale_colour_manual(values = cols)
p3 <- p3 + scale_fill_manual(values = cols)
p3

i4 <- p3 + coord_cartesian(xlim = c(-80, -40), ylim = c(10.5, 15.5)) + theme(legend.position = "none")
i5 <- p3 + coord_cartesian(xlim = c(55, 62.5), ylim = c(-5, 0)) + theme(legend.position = "none")
i6 <- p3 + coord_cartesian(xlim = c(-85, -77), ylim = c(-64, -57)) + theme(legend.position = "none")

p3 + (i6/i5/i4)


#dapc4 - group assignments
dapc4.res <- as.data.frame(dapc4.aln$posterior)
dapc4.res$pop <- aln.f$strata$pop
dapc4.res$ind <- aln.f$strata$seq_code
dapc4.res$spp <- aln.f$strata$spp
dapc4.res$reg <- aln.f$strata$region

dapc4.res.m <-melt(dapc4.res)
colnames(dapc4.res.m)[5:6] <- c("dapc4.group", "hpd")
dapc4.res.m$reg <- factor(dapc4.res.m$reg,levels(as.factor(dapc4.res.m$reg))[c(1, 10, 8, 7, 5, 6, 4, 9, 3, 2)])
dapc4.res.m$p_i <- paste(dapc4.res.m$pop, dapc4.res.m$ind, sep="_")

p4 <- ggplot(dapc4.res.m, aes(x = p_i, y = hpd, fill = dapc4.group))
p4 <- p4 + geom_bar(stat='identity', alpha = 0.8)
p4 <- p4 + scale_fill_manual(values = cols) 
p4 <- p4 + facet_grid(~reg, scales = "free", space = "free")
p4 <- p4 + theme(panel.background = element_rect("white"), 
                 panel.grid.major = element_line("darkgrey"),
                 strip.text.x = element_blank(), 
                 strip.background = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
                 axis.title.x = element_blank(),
                 legend.position = "none")
p4


(((p1 + ((i3+i2 + plot_layout(widths = c(1, 2)))/i1 + plot_layout(heights = c(2,1)))) + plot_layout(widths = c(2, 1))) / (p2 / p4)) + plot_layout(heights = c(5, 1), guides = "collect")
  

#loadings dapc4
set.seed(4)
contrib14 <- loadingplot(dapc4.aln$var.contr, axis = 1, thres = 0.015, lab.jitter = 1)

set.seed(6)
contrib24 <- loadingplot(dapc4.aln$var.contr, axis = 2, thres = 0.015, lab.jitter = 1)





























