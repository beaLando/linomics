library(tidyverse)
library(ape)
library(geiger)
library(phytools)
library(treeio)
library(ggplot2)
library(ggtree)

# Import trees

## nuclear tree
treN <- ape::read.tree(file = "./output/nuclear_lowDepth/genetic_distances_tree.nwk") #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
plot(treN, cex = 0.5)

## plastome tree
treP <- ape::read.tree(file = "./output/plastome/iqtree/without_84/formatted_consensus.nwk.txt") 
treP <- ape::as.phylo(treP)
plot(treP, cex = 0.5)

## tfl1 tree #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
treF <- ape::read.tree(file = "./output/tfl1/tfl1_port_tree/consensus.rooted.tre")


# Root trees & drop tips

#treNr <- ape::root(treN, outgroup = "L88", resolve.root = TRUE, edgelabel = TRUE)
treNr <- ape::drop.tip(treN, tip = c("L84", "mrkdup_HighDepth_IOW2_8", "mrkdup_HighDepth_12_6", "mrkdup_HighDepth_5_1"))
treNr <- ladderize(treNr)
plot(treNr, use.edge.length=FALSE)
nodelabels()

plot(treP, use.edge.length=FALSE)
nodelabels()
#plot(trePr, use.edge.length = FALSE)
#nodelabels()
trePr <- ape::drop.tip(treP, tip = c("Linum_usitatissimum"))
trePr <- ladderize(trePr)
trePr <- ape::rotate(trePr, 96)
#trePr <- ape::root.phylo(trePr, outgroup = "L88", resolve.root = TRUE, edgelabel = TRUE)
plot(trePr, use.edge.length=FALSE)
nodelabels()

plot(treNr, cex = 0.5)
plot(trePr, cex = 0.5)

# Compare trees

both.tre <- cophylo(compute.brlen(trePr, power = 1), compute.brlen(treNr, power = 1), rotate = TRUE)
plot(both.tre, type = "phylogram", 
     fsize = 0.8, 
     link.type = "curved", link.lwd = 3, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"


mrca(trePr)["L88", "L64"]

p.left <- ggtree(trePr, ladderize = F, col = "gray33") #+ #, layout='roundrect'
  #geom_treescale()
  #geom_nodelab() + 
  #geom_nodepoint() + 
  #geom_tippoint() + 
  #geom_tiplab() + 
  #theme_tree2() #+ 
#scale_x_continuous(labels = abs) +
#geom_hilight(
#mapping=aes(subset = node %in% c(38, 48, 58, 36),
#node = node,
#fill = as.factor(node)
#)
#) +
#labs(fill = "Plastome Phylogeny" )

pld <- p.left$data %>% #View()
  mutate(x = if_else(node == 95, x + 0.0083, x)) %>%
  mutate(x = (x*1000)) #%>%
  #mutate(branch.length = if_else(node == 1 | node == 96, branch.length, branch.length*80)) %>%
  #mutate(branch.length = branch.length*2) %>%
  #mutate(branch = if_else(node == 1 | node == 96, branch, branch*1000)) %>%
  #mutate(branch = branch*2)

p.left$data$x <- pld$x
#p.left$data$branch.length <- pld$branch.length
#p.left$data$branch <- pld$branch

#p.left <- p.left + geom_treescale()

p.left

p.right <- ggtree(treNr, ladderize = F, col = "gray33") #+ #, layout='roundrect'
  #geom_nodelab() + 
  #geom_nodepoint() + 
  #geom_tippoint() + 
  #geom_tiplab() + 
  #theme_tree2() #+ 
#scale_x_continuous(labels = abs) +
#geom_hilight(
#mapping=aes(subset = node %in% c(38, 48, 58, 36),
#node = node,
#fill = as.factor(node)
#)
#) +
#labs(fill = "Plastome Phylogeny" )

p.right <- p.right + geom_treescale()

p.right

d1 <- p.left$data
d2 <- p.right$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.2

pp <- p.left + geom_tree(data = d2, col = "gray33") #+
  #geom_treescale() #+
  #ggnewscale::new_scale_fill() #+
#geom_hilight(
#data = d2, 
#mapping = aes( 
#subset = node %in% c(38, 48, 58),
#node=node,
#fill=as.factor(node))
#) +
#labs(fill = "clades for tree in right" ) 

dd.true <- bind_rows(d1, d2) %>% #View()
  filter(!is.na(label)) %>%
  filter(isTip == "TRUE") %>%
  droplevels #%>% View()

dd.false <- bind_rows(d1, d2) %>% 
  filter(!is.na(label)) %>%
  filter(isTip == "FALSE") %>%
  droplevels #%>% View()

pp + geom_line(aes(x, y, group = label), data = dd.true, lwd = 1, col = make.transparent("gray",0.25)) +
  geom_tippoint() + 
  geom_tippoint(data = d2) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), col = make.transparent("blue",0.15)) +
  scale_size_binned(breaks = c(25, 50, 75))
  

 





