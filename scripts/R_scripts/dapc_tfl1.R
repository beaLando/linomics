library(vcfR)
library(poppr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(igraph)
library(tidyverse)

#Import population info
pop.info <- read.csv("./data/lblus_all_info.csv")

#Import MAFFT alignment trimmed in trimaAL
aln <- fasta2DNAbin("./data/Lb_TFL1/lb_all.mft.tal.fasta", quiet = FALSE, chunkSize = 10, snpOnly = FALSE)
aln <- DNAbin2genind(aln, pop = NULL, exp.char = c("a","t","g","c"), polyThres = 1/100)

ploidy(aln) <- 1

#aln.f <- missingno(aln, "loci")
aln.f <- missingno(aln, "geno", cutoff = 0.5)

#Removing 20 genotypes: L96 1281 bp, L91 1281 bp, L84 1281 bp, L68 1281 bp, L88 1281 bp, L35
#1281 bp, L27 1281 bp, L81 1281 bp, L40 1281 bp, L51 1281 bp, L65 1281 bp, L47 1281 bp, L52 1281
#bp, L9 1281 bp, L39 1281 bp, L3 1281 bp, L60 1281 bp, L92 1281 bp, L63 1281 bp, L49 1281 bp

names.fixed <- data.frame(inds = str_remove(indNames(aln.f), " 1281 bp")) %>%
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
  select(-X) %>% 
  filter(seq_code %in% indNames(aln.f)) %>% #str()
  mutate(purpose = case_when(
    purpose == "wild" ~ "wild",
    purpose == "oil" ~ "oil",
    purpose == "fibre" ~ "fiber",
    purpose == "fiber" ~ "fiber",
    purpose == "unknown" ~ "unkn-interm",
    purpose == "intermediate" ~ "unkn-interm"
  )) %>%
  mutate(
    region_purpose = if_else(spp == "usitatissimum", purpose, region),
    LuTFL1 = if_else(is.na(LuTFL1) == TRUE, "new", LuTFL1)) %>%
  distinct(seq_code, .keep_all = TRUE) %>%
  as.data.frame() #%>% View()

aln.f@pop <- as.factor(pop.info2$region_purpose)
aln.f@strata <- pop.info2

# PRODUCE DAPC

#k = max
set.seed(435)
grp <- find.clusters(aln.f, max.n.clust = 20) #PCs = 8, K = 3

table(pop(aln.f), grp$grp)
table.value(table(pop(aln.f), grp$grp), col.lab = paste("k", 1:3))

set.seed(674)
dapc.aln <- dapc(aln.f, var.contrib = TRUE, scale = FALSE, pop =  grp$grp) #PCs = 8, D = 2

scatter(dapc.aln, mstree = TRUE, cell = 0, pch = 1:10, clabel = 0, cstar = 0, legend = T, posi.leg = "right", posi.da = "topleft", ratio.da = 0.2, cleg = 1, scree.pca = T, posi.pca = "bottomleft", ratio.pca = 0.2)
assignplot(dapc.aln)

dapc.plot <- as.data.frame(dapc.aln$ind.coord)
dapc.plot$Group <- dapc.aln$grp
dapc.plot$pop <- aln.f$strata$pop
dapc.plot$ind <- aln.f$strata$seq_code
dapc.plot$spp <- aln.f$strata$spp
dapc.plot$lat <- aln.f$strata$lat
dapc.plot$purpose <- aln.f$strata$purpose
dapc.plot$reg <- aln.f$strata$region
dapc.plot$reg <- factor(dapc.plot$reg,levels(as.factor(dapc.plot$reg))[c(7, 5, 3, 2, 9, 1, 8, 4, 6)])
dapc.plot$regp <- aln.f$strata$region_purpose
dapc.plot$regp <- factor(dapc.plot$regp,levels(as.factor(dapc.plot$regp))[c(1, 10, 8, 7, 4, 5, 3, 9, 11, 6, 2)])
dapc.plot$hap <- aln.f$strata$LuTFL1
head(dapc.plot)

#cols <- brewer.pal(n = nlevels(as.factor(dapc.res.m$dapc.group)), name = "Spectral")
cols <- c("palevioletred4", "darkred", "brown1", "darkorange1", "darkgoldenrod2", "lightgoldenrod",   
"skyblue3", "steelblue4", "darkseagreen4", "gray31")

p1 <- ggplot(dapc.plot, aes(x = LD1, y = LD2, col = reg, fill = reg, shape = purpose))
p1 <- p1 + geom_point(size = 4)
p1 <- p1 + theme_bw(base_size = 20) + theme(legend.position = "left")
p1 <- p1 + scale_shape_manual(values = c(8, 24, 4, 21, 1, 23, 22, 25, 9, 10, 11)) #9, 10, 11, 24 #8, 24, 4, 21, 1, 23, 22, 25, 9, 10, 11  #3, 7, 10, 11
p1 <- p1 + scale_colour_manual(values = cols)
p1 <- p1 + scale_fill_manual(values = cols)
p1


#dapc - group assignments
dapc.res <- as.data.frame(dapc.aln$posterior)
dapc.res$pop <- aln.f$strata$pop
dapc.res$ind <- aln.f$strata$seq_code
dapc.res$spp <- aln.f$strata$spp
dapc.res$regp <- aln.f$strata$region_purpose
dapc.res$lat <- as.character(aln.f$strata$lat)

dapc.res.m <-melt(dapc.res)
colnames(dapc.res.m)[6:7] <- c("dapc.group", "hpd")
dapc.res.m$regp <- factor(dapc.res.m$regp,levels(as.factor(dapc.res.m$regp))[c(1, 10, 8, 7, 4, 5, 3, 9, 11, 6, 2)])
dapc.res.m$p_i <- paste(dapc.res.m$pop, dapc.res.m$ind, sep="_")
dapc.res.m$lat <- as.numeric(dapc.res.m$lat)

dapc.res.m <- dapc.res.m[order(dapc.res.m$lat),]


p2 <- ggplot(dapc.res.m, aes(x = reorder_within(p_i, hpd, regp), y = hpd, fill = dapc.group)) + scale_x_reordered()
p2 <- p2 + geom_bar(stat='identity', alpha = 0.8)
p2 <- p2 + scale_fill_manual(values = c("white", "gray", "black")) 
p2 <- p2 + facet_grid(~regp, scales = "free", space = "free")
p2 <- p2 + theme(panel.background = element_rect("white"), 
                 panel.grid.major = element_line("darkgrey"),
                 strip.background = element_rect(color="darkgrey", fill="white"),
                 strip.text.x = element_text(size = 10, color = "black"),
                 axis.text.x = element_text(size = 3, color = "black", angle = 45),
                 axis.title.x = element_blank(),
                 legend.position = "none")
p2


# PRODUCE HAPLOTYPE NETWORK
aln.f@pop <- as.factor(pop.info2$LuTFL1)
imsn()



































