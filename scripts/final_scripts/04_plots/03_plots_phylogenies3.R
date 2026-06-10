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
treepl_tree_matrix_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/supermx353_linum_sm0_primed_dated.tre")
treepl_tree_astral_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353_linum_sm0_primed_dated.tre")
treepl_tree_matrix_plastid_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/plastid/supermx353_linum_sm0_primed_dated.tre")

# Fix labels
tipsM <- treepl_tree_matrix_gtrg$tip.label
tipsM <- data.frame(tips_old = tipsM, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                 "L. flavum", "L. macraei", "L. strictum", "L. tenuifolium", "L. tenue", 
                                 "L. viscosum", "L. hirsutum", 
                                 "L. grandiflorum", "L. decumbens", "L. marginale", "L. hologynum", 
                                 "L. bienne SW", "L. bienne SE", "L. usitatissimum", "L. bienne NW",  
                                 "L. narbonense D07", "L. narbonense D14", "L. narbonense D08", "L. narbonense D09", 
                                 "L. perenne",  "L. leonii", "L. lewisii"))


treepl_tree_matrix_gtrg$tip.label[match(tipsM$tips_old, treepl_tree_matrix_gtrg$tip.label)] <- tipsM$tips_new
plot(treepl_tree_matrix_gtrg)

tipsA <- treepl_tree_astral_gtrg$tip.label
tipsA <- data.frame(tips_old = tipsA, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                            "L. tenue", "L. flavum", "L. macraei", "L. strictum", "L. tenuifolium", 
                                            "L. perenne",  "L. leonii", "L. lewisii",
                                            "L. grandiflorum", "L. decumbens",
                                            "L. hologynum", "L. marginale",
                                            "L. bienne SW", "L. bienne SE", "L. usitatissimum", "L. bienne NW", 
                                            "L. narbonense D14", "L. narbonense D07", "L. narbonense D08", "L. narbonense D09", 
                                            "L. viscosum", "L. hirsutum"))

treepl_tree_astral_gtrg$tip.label[match(tipsA$tips_old, treepl_tree_astral_gtrg$tip.label)] <- tipsA$tips_new
plot(treepl_tree_astral_gtrg)

treepl_tree_matrix_plastid_gtrg <- drop.tip(treepl_tree_matrix_plastid_gtrg, c("ERR2040369_plastid", "ERR2040381_plastid"))
tipsP <- treepl_tree_matrix_plastid_gtrg$tip.label
tipsP <- data.frame(tips_old = tipsP, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                 "L. flavum", "L. tenue", "L. macraei", "L. strictum", "L. tenuifolium", 
                                 "L. viscosum", "L. hirsutum", 
                                 "L. grandiflorum", "L. hologynum", 
                                 "L. bienne SW", "L. usitatissimum", "L. bienne NW", "L. bienne SE", 
                                 "L. narbonense D07", "L. narbonense D09", "L. narbonense D14", 
                                 "L. perenne",  "L. leonii", "L. lewisii"))

treepl_tree_matrix_plastid_gtrg$tip.label[match(tipsP$tips_old, treepl_tree_matrix_plastid_gtrg$tip.label)] <- tipsP$tips_new
plot(treepl_tree_matrix_plastid_gtrg)

treepl_tree_matrix_gtrg <- drop.tip(treepl_tree_matrix_gtrg, c("E. mesembryanthemifolia", "Phyllanthus spp."))
treepl_tree_astral_gtrg <- drop.tip(treepl_tree_astral_gtrg, c("E. mesembryanthemifolia", "Phyllanthus spp."))
treepl_tree_matrix_plastid_gtrg <- drop.tip(treepl_tree_matrix_plastid_gtrg, c("E. mesembryanthemifolia", "Phyllanthus spp."))

# Plot
p.treepl.m1 <- ggtree(treepl_tree_matrix_gtrg) +
  coord_geo(xlim = c(-40, 4), 
            ylim = c(0, Ntip(treepl_tree_matrix_gtrg)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  theme_tree2()

revts(p.treepl.m1)


p.treepl.a1 <- ggtree(treepl_tree_astral_gtrg) +
  coord_geo(xlim = c(-40, 4), 
            ylim = c(0, Ntip(treepl_tree_astral_gtrg)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  theme_tree2()

revts(p.treepl.a1)

p.treepl.p1 <- ggtree(treepl_tree_matrix_plastid_gtrg) +
  coord_geo(xlim = c(-40, 4), 
            ylim = c(0, Ntip(treepl_tree_matrix_plastid_gtrg)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  theme_tree2()

revts(p.treepl.p1)

revts(p.treepl.p1) + revts(p.treepl.m1) + revts(p.treepl.a1) +
  plot_layout(widths = c(1, 1, 1))


################################ IQTREE - Linum supermatrix & Astral phylogenies ################################
# Import trees
tree_matrix_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/supermx353.phyxed.treefile")
tree_astral_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.phyxed.tre")
tree_matrix_plastid_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/plastid/supermx353.phyxed.treefile")
#tree_matrix_plastid25_gtrg <- ape::read.tree("./data/Angiosperm353/iqtrees/plastid25/supermx353.phyxed.treefile")

# Tree scales
ggtree(tree_matrix_gtrg) +
  geom_tiplab(size = 3) +
  theme_tree2()

ggtree(tree_astral_gtrg) +
  geom_tiplab(size = 3) +
  theme_tree2()

ggtree(tree_matrix_plastid_gtrg) +
  geom_tiplab(size = 3) +
  theme_tree2()

# GTR+G
# Extract coords
left <- ladderize(tree_matrix_gtrg, right = T)
plot(left)
tipsL <- left$tip.label
tipsL <- data.frame(tips_old = tipsL, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                 "L. flavum", "L. macraei", "L. strictum", "L. tenuifolium", "L. tenue", 
                                 "L. viscosum", "L. hirsutum", 
                                 "L. grandiflorum", "L. decumbens", 
                                 "L. marginale", "L. hologynum", 
                                 "L. bienne SW", "L. bienne SE", "L. usitatissimum", "L. bienne NW", 
                                 "L. narbonense D07", "L. narbonense D14", "L. narbonense D08", "L. narbonense D09", 
                                 "L. perenne", "L. leonii", "L. lewisii"))

left$tip.label[match(tipsL$tips_old, left$tip.label)] <- tipsL$tips_new
plot(left)

right <- ladderize(tree_astral_gtrg, right = T)
tipsR <- right$tip.label
tipsR <- data.frame(tips_old = tipsR, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                 "L. tenue", "L. flavum", "L. macraei", "L. strictum", "L. tenuifolium", 
                                 "L. perenne", "L. leonii", "L. lewisii",
                                 "L. grandiflorum", "L. decumbens", 
                                 "L. hologynum", "L. marginale", 
                                 "L. bienne SW", "L. bienne SE", "L. usitatissimum", "L. bienne NW", 
                                 "L. narbonense D14", "L. narbonense D07", "L. narbonense D09", "L. narbonense D08", 
                                 "L. viscosum", "L. hirsutum"))

right$tip.label[match(tipsR$tips_old, right$tip.label)] <- tipsR$tips_new
plot(right)

left2 <- drop.tip(tree_matrix_plastid_gtrg, c("ERR2040369_plastid", "ERR2040381_plastid"))
tips2 <- left2$tip.label
tips2 <- data.frame(tips_old = tips2, 
                    tips_new = c("E. mesembryanthemifolia", "Phyllanthus spp.", 
                                 "L. flavum", "L. tenue", "L. macraei", "L. strictum", "L. tenuifolium", 
                                 "L. viscosum", "L. hirsutum", 
                                 "L. grandiflorum", "L. marginale", 
                                 "L. bienne SW", "L. usitatissimum", "L. bienne NW", "L. bienne SE", 
                                 "L. narbonense D07", "L. narbonense D09", "L. narbonense D14",
                                 "L. perenne",  "L. leonii", "L. lewisii"))

left2$tip.label[match(tips2$tips_old, left2$tip.label)] <- tips2$tips_new
left2 <- ladderize(left2, right = T)
plot(left2)

both.tre1 <- cophylo(left, right, rotate = TRUE)
dev.off()

both.tre2 <- cophylo(left, left2, rotate = TRUE)
dev.off()

plot(both.tre1, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

plot(both.tre2, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()

left <- both.tre1[[1]][[1]]
right <- both.tre1[[1]][[2]]
left2 <- both.tre2[[1]][[2]]

# Plot trees face to face in ggtree
p.left <- ggtree(left, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
#geom_treescale(width =  0.1) +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() + #95 is blue+lightblue # 98 is yellow+red #100 is yellow #95 is red #78 is light blue
# scale_x_continuous(labels = abs)

p.right <- ggtree(right, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
#geom_treescale(width =  0.1) +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() +
# scale_x_continuous(labels = abs)

p.left2 <- ggtree(left2, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
#geom_treescale(width = 0.01) +
# geom_nodelab() +
# geom_nodepoint() +
# geom_tippoint() +
# geom_tiplab() +
# theme_tree2() +
# scale_x_continuous(labels = abs)

pld <- p.left$data %>% #View()
  mutate(x = (x*28))

p.left$data$x <- pld$x

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

p.right <- p.right + geom_treescale()

p.right + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

pld2 <- p.left2$data %>% #View()
  mutate(x = (x*55)) %>%
  mutate(y = (y*1.15))

p.left2$data$x <- pld2$x
p.left2$data$y <- pld2$y

p.left2 + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

d1 <- p.left$data
d2 <- p.right$data
d3 <- p.left2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 5
d3$x <- -(max(d3$x) - d3$x + min(d1$x))

pp <- p.left + geom_tree(data = d2, col = "gray33") + geom_tree(data = d3, col = "gray33")

dd.true <- bind_rows(d1, d2, d3) %>% #View()
  filter(!is.na(label)) %>%
  filter(isTip == "TRUE") %>%
  droplevels #%>% View()

dd.false <- bind_rows(d1, d2, d3) %>% 
  filter(isTip == "FALSE") %>% 
  mutate(label = if_else(as.numeric(label) <= 1, as.numeric(label)*100, as.numeric(label))) %>%
  filter(!is.na(label)) %>%
  droplevels #%>% View()

pp + geom_line(aes(x, y, group = label), data = dd.true, lwd = 1, col = make.transparent("gray",0.25)) +
  geom_tippoint() + 
  geom_tippoint(data = d2) + 
  geom_tippoint(data = d3) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d3, hjust = -0.3, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
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
treP <- ape::read.tree(file = "./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/all_plastomes_gtrg.aligned.trimmed.phyxed.treefile")
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
  mutate(x = (x*2000)) #2000

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
beast_tree_matrix <- read.beast("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_superMX5_10pc_mcc.tre")
beast_tree_astral <- read.beast("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_Astral4_mcc.tre")
beast_tree_plastid <- read.beast("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_plastid_10pc_MCC.tre")

beast_tree_matrix@data$height_median <- ifelse(beast_tree_matrix@data$height_median < 1, NA, beast_tree_matrix@data$height_median)
beast_tree_astral@data$height_median <- ifelse(beast_tree_astral@data$height_median < 0.1, NA, beast_tree_astral@data$height_median)
beast_tree_plastid@data$height_median <- ifelse(beast_tree_plastid@data$height_median < 1, NA, beast_tree_plastid@data$height_median)

# Plot including Narbonense
p.beast.m <- revts(ggtree(beast_tree_matrix)) +
  coord_geo(xlim = c(-18, 2),
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

p.beast.a <- revts(ggtree(beast_tree_astral)) +
  coord_geo(xlim = c(-2.5, 0.3),
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

p.beast.p <- revts(ggtree(beast_tree_plastid)) +
  coord_geo(xlim = c(-20, 2),
            ylim = c(-2, Ntip(beast_tree_plastid)),
            dat = list("epochs", "periods"),
            neg = TRUE, abbrv = FALSE, center_end_labels = TRUE,
            fill = "transparent") +
  geom_tiplab(size = 1.5) + 
  geom_range(range = "height_0.95_HPD", 
             center = "height_median",
             col = "blue",
             alpha = 0.5) +
  theme_tree2()


p.beast.p + p.beast.m + p.beast.a


# Also plot parameters
## Supermx calibration
beast_tree_matrix <- load.trees("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_superMX5.(time).trees", 
                                logfile = "./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_superMX5.log", 
                                format = "beast")
matrix.rwty <- makeplot.all.params(beast_tree_matrix, burnin=10000)

### Save histograms (2:29 are traces)
p.list <- list()

for (i in 1:30){
  p.list[[i]] <- matrix.rwty[[i]]$density.plot + 
    scale_fill_manual(values = c("gray","black")) + 
    labs(title = NULL) +
    theme_classic() + 
    theme(legend.position = "none")
}

p.list <- p.list[-c(8,9,28)]

grid.arrange(grobs = p.list, ncol = 9)
remove(beast_tree_matrix)
remove(matrix.rwty)

## Astral calibration
beast_tree_astral <- load.trees("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_Astral4.(time).trees", 
                                logfile = "./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_Astral4.log", 
                                format = "beast")
astral.rwty <- makeplot.all.params(beast_tree_astral, burnin=10000)

### Save histograms (2:29 are traces)
p.list <- list()

for (i in 1:30){
  p.list[[i]] <- astral.rwty[[i]]$density.plot + 
    scale_fill_manual(values = c("gray","black")) + 
    labs(title = NULL) +
    theme_classic() + 
    theme(legend.position = "none")
}

p.list <- p.list[-c(8,9,28)]

grid.arrange(grobs = p.list, ncol = 9)
remove(beast_tree_astral)
remove(astral.rwty)

## Plastid calibration
beast_tree_plastid <- load.trees("./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_plastid.(time).trees", 
                                logfile = "./data/Lb_plastomes_raw/Consensus/aln/all/newRes2026/run_plastid.log", 
                                format = "beast")
plastid.rwty <- makeplot.all.params(beast_tree_plastid, burnin=10000)

### Save histograms (2:29 are traces)
p.list <- list()

for (i in 1:30){
  p.list[[i]] <- plastid.rwty[[i]]$density.plot + 
    scale_fill_manual(values = c("gray","black")) + 
    labs(title = NULL) +
    theme_classic() + 
    theme(legend.position = "none")
}

p.list <- p.list[-c(8,9,28)]

grid.arrange(grobs = p.list, ncol = 9)
remove(beast_tree_plastid)
remove(plastid.rwty)
