library(tidyverse)
library(ape)
library(geiger)
library(phytools)
library(treeio)
library(ggplot2)
library(ggtree)

# Import trees

## nuclear tree
treAst <- ape::read.tree(file = "./output/Angiosperm353/Astral/gtr_i_g/astral353_rootedITOL.nwk") #I rooted tree in itol because somehow iqtree output does not work with ggtree otherwise
treAst #check if rooted, else: <- root(treAst, outgroup = "Euphorbia_mesembryanthemifolia", resolve.root = TRUE)
plot(treAst, cex = 0.5)

## plastome tree
treMtx <- ape::read.tree(file = "./output/Angiosperm353/iqtree/gtr_i_g/partitions_rooted/supermx_rootedITOL.nwk")
treMtx #check if rooted, else: <- root(treMtx, outgroup = "Euphorbia_mesembryanthemifolia", resolve.root = TRUE)
plot(treMtx, cex = 0.5, show.node.label = TRUE)



left <- ladderize(treMtx, right = T)
plot(left)
add.scale.bar()

right <- ladderize(treAst, right = T)
plot(right)
add.scale.bar()

both.tre <- cophylo(left, right, rotate = TRUE)



dev.off()

h1<-max(nodeHeights(left))
h2<-max(nodeHeights(right))

plot(both.tre, type = "phylogram", 
     fsize = 0.5,
     link.type = "curved", link.lwd = 1.5, link.lty = "solid",
     link.col = make.transparent("blue",0.25)) #, link.type = "curved"

dev.off()

left <- read.nhx(file = "./output/Angiosperm353/iqtree/gtr_i_g/partitions_rooted/supermx_rootedITOL.nwk")
right <- read.nhx(file = "./output/Angiosperm353/Astral/gtr_i_g/astral353_rootedITOL.nwk") 


# Re-plot in ggplot

p.left <- ggtree(left, ladderize = TRUE, col = "gray33") #+ geom_text(aes(label=node))#+ #geom_treescale()#, layout='roundrect'
p.left <- p.left %>% flip(9,10) %>% flip(3,4)
# geom_nodelab() + 
# geom_nodepoint() + 
# geom_tippoint() + 
# geom_tiplab() + 
# theme_tree2() #+ 
#scale_x_continuous(labels = abs) +
#geom_hilight(
#mapping=aes(subset = node %in% c(38, 48, 58, 36),
#node = node,
#fill = as.factor(node)
#)
#) +
#labs(fill = "Plastome Phylogeny" )

pld <- p.left$data %>% #View()
  #mutate(x = if_else(node == 95, x + 0.0083, x)) %>%
  mutate(x = (x*(100/3))) #%>%
#mutate(branch.length = if_else(node == 1 | node == 96, branch.length, branch.length*80)) %>%
#mutate(branch.length = branch.length*2) %>%
#mutate(branch = if_else(node == 1 | node == 96, branch, branch*1000)) %>%
#mutate(branch = branch*2)

p.left$data$x <- pld$x
#p.left$data$branch.length <- pld$branch.length
#p.left$data$branch <- pld$branch

#p.left <- p.left + geom_treescale()

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

p.right <- ggtree(right, ladderize = TRUE, col = "gray33") #+ geom_treescale() #, layout='roundrect'
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

p.right + geom_text(aes(label=node), hjust=-.3) + geom_treescale()


d1 <- p.left$data
d2 <- p.right$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 6

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
  geom_nodepoint(aes(size = B), col = "gray", alpha = 0.5) +
  #geom_tippoint() + 
  #geom_tippoint(data = d2) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), col = make.transparent("blue",0.15)) +
  scale_size_binned(breaks = c(25, 50, 75, 100)) + geom_treescale(x = 0, y = 25)



pp + geom_line(aes(x, y, group = label), data = dd.true, lwd = 1, col = make.transparent("gray",0.25)) +
  geom_nodepoint(aes(size = B), col = "gray", alpha = 0.5) +
  #geom_tippoint() + 
  #geom_tippoint(data = d2) + 
  geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), col = make.transparent("blue",0.15)) +
  scale_size_binned(breaks = c(25, 50, 75)) + geom_treescale(x = 0, y = 25) +
  coord_cartesian(x = c(11.5, 19), y = c(23, 26))




