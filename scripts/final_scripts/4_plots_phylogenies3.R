################################### LIBRARIES ###################################
library(tidyverse)
library(ape)
library(geiger)
library(phytools)
library(treeio)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ggtree)
library(deeptime)
library(rwty)
library(raster)


# ################################### DENSITREE - Linum supermatrix phylogeny ###################################
# # Import trees
# ## nuclear tree
# btrees <- treeio::read.tree("./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.ufboot")
# 
# ggdensitree(btrees, 
#             alpha = .3, colour = 'steelblue') + 
#   geom_tiplab(size = 3) + 
#   hexpand(.35)


################################ TREEPL - Linum supermatrix & Astral phylogenies dated ################################
# Import trees
treepl_tree_matrix <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353_linum_sm0_primed_dated.tre")
treepl_tree_astral <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_linum_sm0_primed_dated.tre")

# Plot including Narbonense
p.treepl.m <- ggtree(treepl_tree_matrix) +
  coord_geo(xlim = c(-90, 4), 
            ylim = c(0, Ntip(treepl_tree_matrix)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  theme_tree2()

revts(p.treepl.m)


p.treepl.a <- ggtree(treepl_tree_astral) +
  coord_geo(xlim = c(-50, 4), 
            ylim = c(0, Ntip(treepl_tree_astral)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  theme_tree2()

revts(p.treepl.a)

revts(p.treepl.m) + revts(p.treepl.a) +
  plot_layout(widths = c(1, 0.6))


################################ IQTREE - Linum supermatrix & Astral phylogenies ################################
# Import trees
treepl_tree_matrix <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.phyxed.contree")
treepl_tree_astral <- ape::read.tree("./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre")

# Tree scales
ggtree(treepl_tree_matrix) +
  geom_tiplab(size = 3) +
  theme_tree2()

ggtree(treepl_tree_astral) +
  geom_tiplab(size = 3) +
  theme_tree2()

# Extract coords
left <- ladderize(treepl_tree_matrix, right = T)
plot(left)

right <- ladderize(treepl_tree_astral, right = T)
plot(right)

both.tre <- cophylo(left, right, rotate = TRUE)
dev.off()

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()

left <- both.tre[[1]][[1]]
right <- both.tre[[1]][[2]]

# Plot trees face to face in ggtree
p.left <- ggtree(left, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
# geom_treescale() +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() + #95 is blue+lightblue # 98 is yellow+red #100 is yellow #95 is red #78 is light blue
# scale_x_continuous(labels = abs)

p.right <- ggtree(right, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
# geom_treescale() +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() +
# scale_x_continuous(labels = abs)

pld <- p.left$data %>% #View()
  mutate(x = (x*30))

p.left$data$x <- pld$x

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

p.right <- p.right + geom_treescale()

p.right + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

d1 <- p.left$data
d2 <- p.right$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 5

pp <- p.left + geom_tree(data = d2, col = "gray33")

dd.true <- bind_rows(d1, d2) %>% #View()
  filter(!is.na(label)) %>%
  filter(isTip == "TRUE") %>%
  droplevels #%>% View()

dd.false <- bind_rows(d1, d2) %>% 
  filter(isTip == "FALSE") %>% 
  mutate(label = if_else(as.numeric(label) <= 1, as.numeric(label)*100, as.numeric(label))) %>%
  filter(!is.na(label)) %>%
  droplevels #%>% View()

pp + geom_line(aes(x, y, group = label), data = dd.true, lwd = 1, col = make.transparent("gray",0.25)) +
  geom_tippoint() + 
  geom_tippoint(data = d2) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), 
             col = make.transparent("blue",0.15)) +
  scale_size_binned(breaks = c(30, 60, 90)) +
  geom_treescale(x = 1, y = 1)


################################### IQTREE - Chloroplast genes phylogeny w/ L. villarianum ###################################
# Import trees
## nuclear tree
treCgs <- ape::read.tree(file = "./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx.phyxed.treefile")
plot(treCgs, cex = 0.5)

treCgs <- ladderize(treCgs, right = T)
plot(treCgs)

both.tre <- cophylo(treCgs, treCgs, rotate = TRUE)
dev.off()

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()

treCgs <- both.tre[[1]][[1]]

ggtree(treCgs, ladderize = FALSE, col = "gray33") + geom_tiplab(hjust = -0.25, size = 2.5) +
  geom_tippoint() + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_nodelab(hjust = 1.5, size = 3, col = "blue") #+
#scale_size_binned(breaks = c(30, 60, 90)) +
#geom_treescale(x = 1, y = 1)


################################### IQTREE - Chloroplast and Nuclear Lb phylogenies ###################################
# Import trees
## nuclear tree
treN <- ape::read.tree(file = "./data/Lb_plastomes_raw/Consensus/genetic_distances_tree_nuclear.phyxed.nwk") #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
plot(treN, cex = 0.5)

## plastome tree
treP <- ape::read.tree(file = "./data/Lb_plastomes_raw/Consensus/aln/all/all_plastomes.aligned.trimmed.fasta.phyxed.contree")
plot(treP, cex = 0.5)

left <- ladderize(treP, right = T)
plot(left)
left2 <- drop.tip(left, "L88")
plot(left2)

right <- ladderize(treN, right = T)
plot(right)
right2 <- drop.tip(right, c("L88", "L84"))
plot(right2)

both.tre <- cophylo(left2, right2, rotate = TRUE)
dev.off()

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()

left <- both.tre[[1]][[1]]
right <- both.tre[[1]][[2]]


# Plot trees face to face in ggtree

p.left <- ggtree(left, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
# geom_treescale() +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() + #95 is blue+lightblue # 98 is yellow+red #100 is yellow #95 is red #78 is light blue
# scale_x_continuous(labels = abs)

p.right <- ggtree(right, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
# geom_treescale() +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() +
# scale_x_continuous(labels = abs)

pld <- p.left$data %>% #View()
  mutate(x = (x*2000))

p.left$data$x <- pld$x

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

p.right <- p.right + geom_treescale()

p.right + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

d1 <- p.left$data
d2 <- p.right$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.2

pp <- p.left + geom_tree(data = d2, col = "gray33")

dd.true <- bind_rows(d1, d2) %>% #View()
  filter(!is.na(label)) %>%
  filter(isTip == "TRUE") %>%
  droplevels #%>% View()

dd.false <- bind_rows(d1, d2) %>% 
  filter(isTip == "FALSE") %>% 
  mutate(label = if_else(as.numeric(label) <= 1, as.numeric(label)*100, as.numeric(label))) %>%
  filter(!is.na(label)) %>%
  droplevels #%>% View()

pp1 <- pp + geom_line(aes(x, y, group = label), data = dd.true, lwd = 1, col = make.transparent("gray",0.25)) +
  geom_tippoint() + 
  geom_tippoint(data = d2) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), 
             col = make.transparent("blue",0.15)) +
  scale_size_binned(breaks = c(30, 60, 90)) +
  geom_treescale(x = 1, y = 1)


# Plot trees + niche model map
# raster.present <- raster("./data/niche_model/maxent_preds_lb.tif")
# raster.mask <- raster("./data/niche_model/maxent_preds_lb_thr.tif")
# 
# plot(raster.present)
# plot(raster.mask, col = "darkgray")
# 
# raster_spdf <- as(raster.present, "SpatialPixelsDataFrame")
# raster_df <- as.data.frame(raster_spdf)
# colnames(raster_df) <- c("value", "x", "y")
# 
# if (require("maps")) {
#   ggplot() +
#     geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) + 
#     borders("world", xlim = c(-25, 70), ylim = c(20, 80), col = "black") +
#     coord_cartesian(xlim = c(-25, 63), ylim = c(25, 60)) +
#     scale_fill_gradient(low = "white", high = "black",
#                         limits = c(0, 1)) +
#     theme_minimal()
# }

## Exrtact labels chloroplast
plot(left2, cex = 0.5)
nodelabels(cex = 0.5)

chl.east <- left2$tip.label[phangorn::Descendants(left2, node = 145, type = "tips")[[1]]]
chl.nw <- left2$tip.label[phangorn::Descendants(left2, node = 150, type = "tips")[[1]]]
chl.cnt <- left2$tip.label[phangorn::Descendants(left2, node = 134, type = "tips")[[1]]]
chl.sw <- left2$tip.label[phangorn::Descendants(left2, node = 97, type = "tips")[[1]]]

chl.east <- setdiff(chl.east, chl.nw)

assigned.chl <- data.frame(
  chloroplast_lineage = c(rep("south-east", length(chl.east)),
                          rep("north-west", length(chl.nw)),
                          rep("central", length(chl.cnt)),
                          rep("south-west", length(chl.sw))),
  accession = c(chl.east, chl.nw, chl.cnt, chl.sw))

## Extract labels nuclear
plot(right2, cex = 0.5)
nodelabels(cex = 0.5)

nuc.east <- right2$tip.label[phangorn::Descendants(right2, node = 97, type = "tips")[[1]]]
nuc.nw <- right2$tip.label[phangorn::Descendants(right2, node = 109, type = "tips")[[1]]]
nuc.cnt <- right2$tip.label[phangorn::Descendants(right2, node = 183, type = "tips")[[1]]]
nuc.sw <- right2$tip.label[phangorn::Descendants(right2, node = 153, type = "tips")[[1]]]

nuc.east <- setdiff(nuc.east, nuc.nw)

assigned.nuc <- data.frame(
  nuclear_lineage = c(rep("south-east", length(nuc.east)),
                      rep("north-west", length(nuc.nw)),
                      rep("central", length(nuc.cnt)),
                      rep("south-west", length(nuc.sw))),
  accession = c(nuc.east, nuc.nw, nuc.cnt, nuc.sw))

## Merge labels with accessions info & plot with raster
assigned.lineages <- assigned.chl %>%
  full_join(.,
            assigned.nuc,
            by = "accession") %>%
  dplyr::select(c(accession, chloroplast_lineage, nuclear_lineage)) %>%
  as.data.frame()

pop.info <- xlsx::read.xlsx("./data/populations_locations_1_coordFixed.xlsx", sheetName = "pop_info_plastome")
str(pop.info)

pop.info <- pop.info %>%
  left_join(.,
            assigned.lineages,
            join_by("seq_code"=="accession")) %>%
  as.data.frame() #%>% View()

pop.info <- pop.info %>% 
  filter(!(seq_code %in% c("L88", "L84", "L16"))) %>% 
  filter(!is.na(chloroplast_lineage)) %>%
  filter(!is.na(nuclear_lineage)) %>%
  droplevels()

if (require("maps")) {
  ggplot() +
    #geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) + 
    borders("world", xlim = c(-25, 70), ylim = c(20, 80), col = "black") +
    geom_point(data = pop.info, aes(x = as.numeric(lon), y = as.numeric(lat), 
                                    shape = spp, 
                                    col = chloroplast_lineage), 
               size = 6) +
    scale_colour_manual(values = c("darkgoldenrod1", "lightskyblue", "dodgerblue3", "red")) +
    ggnewscale::new_scale_color() +
    geom_point(data = pop.info, aes(x = as.numeric(lon), y = as.numeric(lat), 
                                    shape = spp, 
                                    col = nuclear_lineage), 
               size = 4) +
    scale_colour_manual(values = c("darkgoldenrod1", "lightskyblue", "dodgerblue3", "red")) +
    coord_cartesian(xlim = c(-25, 67), ylim = c(25, 60)) +
    scale_fill_gradient(low = "white", high = "black",
                        limits = c(0, 1)) +
    labs(x = "Longitude", y = "Latitude", fill = "Suitability", shape = "Species") +
    theme_classic(base_size = 18) +
    theme(legend.position = "none")
}

View(pop.info)


################################### BEAST - Chloroplast Lb phylogeny ###################################
# Import trees
beast_tree_matrix <- read.beast("./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_superMX_treeAnn.(time).trees")
beast_tree_astral <- read.beast("./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_Astral_treeAnn.(time).trees")

beast_tree_matrix@data$height_median <- ifelse(beast_tree_matrix@data$height_median < 1, NA, beast_tree_matrix@data$height_median)
beast_tree_astral@data$height_median <- ifelse(beast_tree_astral@data$height_median < 1, NA, beast_tree_astral@data$height_median)

# Plot including Narbonense
p.beast.m <- ggtree(beast_tree_matrix) +
  coord_geo(xlim = c(-21, 2), 
            ylim = c(-2, Ntip(beast_tree_matrix)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  geom_range(range = "height_0.95_HPD", 
             center = "height_median",
             col = "blue",
             alpha = 0.5) +
  theme_tree2()

revts(p.beast.m)


p.beast.a <- ggtree(beast_tree_astral) +
  coord_geo(xlim = c(-21, 2), 
            ylim = c(-2, Ntip(beast_tree_astral)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  geom_range(range = "height_0.95_HPD", 
             center = "height_median",
             col = "blue",
             alpha = 0.5) +
  theme_tree2()

revts(p.beast.a)


# Plot without narbonense
beast_tree_matrix_d <- drop.tip(beast_tree_matrix, "L88")
beast_tree_astral_d <- drop.tip(beast_tree_astral, "L88")

p.beast.md <- ggtree(beast_tree_matrix_d) +
  coord_geo(xlim = c(-7.5, 0.7), 
            ylim = c(-2, Ntip(beast_tree_matrix_d)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  geom_range(range = "height_0.95_HPD", 
             center = "height_median",
             col = "blue",
             alpha = 0.5) +
  theme_tree2()

revts(p.beast.md)


p.beast.ad <- ggtree(beast_tree_astral_d) +
  coord_geo(xlim = c(-7.5, 0.7), 
            ylim = c(-2, Ntip(beast_tree_astral_d)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  geom_range(range = "height_0.95_HPD", 
             center = "height_median",
             col = "blue",
             alpha = 0.5) +
  theme_tree2()

revts(p.beast.ad)


# Also plot parameters
## Supermx calibration
beast_tree_matrix <- load.trees("./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_superMX.(time).trees.txt", logfile = "./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_superMX.log.txt", format = "beast")
matrix.rwty <- analyze.rwty(beast_tree_matrix, burnin = 1000)

### Save histograms (2:29 are traces)
p.list <- list()

for (i in 1:28){
  j <- c(2:29)[i]
  p.list[[i]] <- matrix.rwty[[j]]$density.plot
}

grid.arrange(grobs = p.list, ncol = 7)
remove(beast_tree_matrix)
remove(matrix.rwty)

## Astral calibration
beast_tree_astral <- load.trees("./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_Astral.(time).trees.txt", logfile = "./data/Lb_plastomes_raw/Consensus/aln/all_plastomes.aligned.trimmed_Astral.log.txt", format = "beast")
astral.rwty <- analyze.rwty(beast_tree_astral, burnin = 1000)

### Save histograms (2:29 are traces)
p.list <- list()

for (i in 1:28){
  j <- c(2:29)[i]
  p.list[[i]] <- astral.rwty[[j]]$density.plot
}

grid.arrange(grobs = p.list, ncol = 7)
remove(beast_tree_astral)
remove(astral.rwty)