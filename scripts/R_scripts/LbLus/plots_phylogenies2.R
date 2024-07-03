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




################################################ PLOT WITH TIPS ########################################################

p.left <- ggtree(left, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
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
  #mutate(x = if_else(node == 95, x + 0.0083, x)) %>%
  mutate(x = (x*2000)) #%>%
#mutate(branch.length = if_else(node == 1 | node == 96, branch.length, branch.length*80)) %>%
#mutate(branch.length = branch.length*2) %>%
#mutate(branch = if_else(node == 1 | node == 96, branch, branch*1000)) %>%
#mutate(branch = branch*2)

p.left$data$x <- pld$x
#p.left$data$branch.length <- pld$branch.length
#p.left$data$branch <- pld$branch

#p.left <- p.left + geom_treescale()

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

collapse(p.left, 151, "max", color = "red", fill = "red") %>%
  collapse(140, "max", color = "goldenrod", fill = "goldenrod") %>%
  collapse(135, "max", color = "darkblue", fill = "darkblue") %>%
  collapse(97, "max", color = "lightblue", fill = "lightblue")

p.right <- ggtree(right, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
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

p.right2 <- groupClade(right, c(96, 153, 97, 110))
ggtree(p.right2, ladderize = FALSE, aes(col = group)) +
  #geom_point(aes(size = as.numeric(label)), col = make.transparent("blue",0.15)) +
  scale_colour_manual(values = c("red", "gold", "dodgerblue4", 'deepskyblue1')) +
  geom_treescale()

p.right2 + 

collapse(p.right, 160, "max", color = "red", fill = "red") %>%
  collapse(153, "max", color = "goldenrod", fill = "goldenrod") %>%
  collapse(148, "max", color = "darkblue", fill = "darkblue") %>%
  collapse(99, "max", color = "lightblue", fill = "lightblue")

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
  scale_size_binned(breaks = c(25, 50, 75)) + geom_treescale(x = 0, y = 1)





################################################ PLOT WITH COLLAPSED LINEAGES ########################################################

p.left <- ggtree(left, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'#+
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
  #mutate(x = if_else(node == 95, x + 0.0083, x)) %>%
  mutate(x = (x*2000)) #%>%
#mutate(branch.length = if_else(node == 1 | node == 96, branch.length, branch.length*80)) %>%
#mutate(branch.length = branch.length*2) %>%
#mutate(branch = if_else(node == 1 | node == 96, branch, branch*1000)) %>%
#mutate(branch = branch*2)

p.left$data$x <- pld$x
#p.left$data$branch.length <- pld$branch.length
#p.left$data$branch <- pld$branch

#p.left <- p.left + geom_treescale()

p.left + geom_text(aes(label=node), hjust=-.3) + geom_treescale()

collapse(p.left, 151, "max", color = "red", fill = "red") %>%
  collapse(140, "max", color = "goldenrod", fill = "goldenrod") %>%
  collapse(135, "max", color = "darkblue", fill = "darkblue") %>%
  collapse(97, "max", color = "lightblue", fill = "lightblue")

p.right <- ggtree(right, ladderize = FALSE, col = "gray33") #+ #, layout='roundrect'
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

collapse(p.right, 160, "max", color = "red", fill = "red") %>%
  collapse(153, "max", color = "goldenrod", fill = "goldenrod") %>%
  collapse(148, "max", color = "darkblue", fill = "darkblue") %>%
  collapse(99, "max", color = "lightblue", fill = "lightblue")

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
  #geom_tiplab(hjust = -0.25, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  #geom_tiplab(data = d2, hjust = 1.7, size = 2.5) + #geom = "label", fill = "white", size = 2, label.padding = unit(0.1, "lines")
  geom_point(data = dd.false, aes(size = as.numeric(label)), col = make.transparent("blue",0.15)) +
  geom_hilight(mapping = aes( 
    subset = node %in% c(151, 140, 135, 97),
    node=node,
    fill=as.factor(node))) +
  geom_hilight(data = d2, mapping = aes( 
    subset = node %in% c(160, 153, 148, 99),
    node=node,
    fill=as.factor(node))) +
  scale_size_binned(breaks = c(25, 50, 75)) + geom_treescale(x = 0, y = 1)












