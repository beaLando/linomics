# LIBS

library(readxl)
library(tidyverse)

library(ggplot2)
library(ggnewscale)

library(raster)
library(maps)
library(rgdal)
library(ecospat)

library(heatmaply)
library(HH)
library(factoextra)

library(phytools)
library(TreeTools)


# DATA

## Coordinates
dat <- read_excel("G:/Shared drives/PHDthesis/5_Phylogeography/v2/Figures/Table2_Lbienne_Lusitatissimum_samples.xlsx") 
glimpse(dat)

## Plastome Phylogeny Sample Assignation 
dat <- dat %>% 
  mutate(.,
         Plastome = recode(seq_code,
                           "L40" = "South West",
                           "L39" = "South West",
                           "L35" = "South West",
                           "L34" = "South West",
                           "L37" = "South West",
                           "L36" = "South West",
                           "L09" = "South West",
                           "L32" = "South West",
                           "L15" = "South West",
                           "L14" = "South West",
                           "L08" = "South West",
                           "L33" = "South West",
                           "L38" = "South West",
                           "L89" = "South West",
                           "L26" = "South West",
                           "L23" = "South West",
                           "L22" = "South West",
                           "L29" = "South West",
                           "L25" = "South West",
                           "L28" = "South West",
                           "L62" = "South West",
                           "L57" = "South West",
                           "L30" = "South West",
                           "L31" = "South West",
                           "L19" = "South West",
                           "L17" = "South West",
                           "L21" = "South West",
                           "L20" = "South West",
                           "L90" = "South West",
                           "L27" = "South West",
                           "L24" = "South West",
                           "L03" = "South West",
                           "L91" = "South West",
                           "L74" = "South West",
                           "L70" = "South West",
                           "L73" = "South West",
                           "L41" = "South West",
                           "L18" = "CM",
                           "L02" = "CM",
                           "L01" = "CM",
                           "L69" = "CM",
                           "L68" = "CM",
                           "L87" = "CM",
                           "L05" = "CM",
                           "L04" = "CM",
                           "L50" = "CM",
                           "L48" = "CM",
                           "L49" = "CM",
                           "L47" = "CM",
                           "L72" = "East",
                           "L71" = "East",
                           "L86" = "East",
                           "L75" = "East",
                           "L81" = "East",
                           "L13" = "East",
                           "L12" = "East",
                           "L92" = "East",
                           "L94" = "East",
                           "L95" = "East",
                           "L96" = "East",
                           "L82" = "East",
                           "L93" = "East",
                           "L80" = "East",
                           "L51" = "North West",
                           "L83" = "North West",
                           "L77" = "North West",
                           "L76" = "North West",
                           "L79" = "North West",
                           "L78" = "North West",
                           "L63" = "North West",
                           "L06" = "North West",
                           "L67" = "North West",
                           "L59" = "North West",
                           "L85" = "North West",
                           "L58" = "North West",
                           "L55" = "North West",
                           "L54" = "North West",
                           "L65" = "North West",
                           "L64" = "North West",
                           "L46" = "North West",
                           "L11" = "North West",
                           "L66" = "North West",
                           "L61" = "North West",
                           "L60" = "North West",
                           "L56" = "North West",
                           "L42" = "North West",
                           "L10" = "North West",
                           "L44" = "North West",
                           "L52" = "East",
                           "L45" = "East",
                           "L53" = "East",
                           "L07" = "East",
                           "L43" = "East",
                           "Lus" = "South West",
                           "L84" = "Outgroups", 
                           "L88" = "Outgroups"))


# Nuclear Phylogeny Sample Assignation

dat <- dat %>% 
  mutate(.,
         Nuclear = recode(seq_code,
                          "L89" = "South West",
                          "L09" = "South West",
                          "L08" = "South West",
                          "L15" = "South West",
                          "L14" = "South West",
                          "L40" = "South West",
                          "L39" = "South West",
                          "L38" = "South West",
                          "L37" = "South West",
                          "L35" = "South West",
                          "L34" = "South West",
                          "L32" = "South West",
                          "L36" = "South West",
                          "L33" = "South West",
                          "L30" = "South West",
                          "L31" = "South West",
                          "L25" = "South West",
                          "L29" = "South West",
                          "L28" = "South West",
                          "L23" = "South West",
                          "L22" = "South West",
                          "L21" = "South West",
                          "L20" = "South West",
                          "L26" = "South West",
                          "L27" = "South West",
                          "L62" = "South West",
                          "L57" = "South West",
                          "L19" = "South West",
                          "L17" = "South West",
                          "L90" = "South West",
                          "L18" = "CM",
                          "L03" = "CM",
                          "L02" = "CM",
                          "L01" = "CM",
                          "L41" = "CM",
                          "L69" = "CM",
                          "L68" = "CM",
                          "L91" = "East",
                          "L74" = "East",
                          "L70" = "East",
                          "L73" = "East",
                          "L72" = "East",
                          "L71" = "East",
                          "L94" = "East",
                          "L93" = "East",
                          "L95" = "East",
                          "L77" = "East",
                          "L76" = "East",
                          "L92" = "East",
                          "L80" = "East",
                          "L79" = "East",
                          "L78" = "East",
                          "L81" = "East",
                          "L82" = "East",
                          "L83" = "East",
                          "L86" = "East",
                          "L75" = "East",
                          "L96" = "East",
                          "L13" = "North West",
                          "L12" = "North West",
                          "L87" = "North West",
                          "L46" = "North West",
                          "L45" = "North West",
                          "L05" = "North West",
                          "L04" = "North West",
                          "L50" = "North West",
                          "L48" = "North West",
                          "L51" = "North West",
                          "L52" = "North West",
                          "L59" = "North West",
                          "L56" = "North West",
                          "L49" = "North West",
                          "L47" = "North West",
                          "L67" = "North West",
                          "L66" = "North West",
                          "L42" = "North West",
                          "L85" = "North West",
                          "L55" = "North West",
                          "L54" = "North West",
                          "L53" = "North West",
                          "L63" = "North West",
                          "L61" = "North West",
                          "L60" = "North West",
                          "L11" = "North West",
                          "L10" = "North West",
                          "L43" = "North West",
                          "L58" = "North West",
                          "L06" = "North West",
                          "L07" = "North West",
                          "L64" = "North West",
                          "L65" = "North West",
                          "L44" = "North West",
                          "Lus" = NA_character_,
                          "L84" = "Outgroups", 
                          "L88" = "Outgroups"))




table(dat$Plastome)
table(dat$Nuclear)

rmvd.nlb <- dat %>%
  filter(spp %in% c("usitatissimum", "narbonense", "NA")) %>%
  as.data.frame()

rmvd.nlb <- rmvd.nlb$seq_code

dat <- dat %>% 
  filter(Plastome != "Outgroups") %>% 
  filter(seq_code != "Lus") %>% 
  mutate_at(.,
            vars(longitude, latitude),
            as.numeric) %>% 
  filter(spp != "usitatissimum")


## Climate

countries <- readOGR("./data/niche_model/countries/world_countries.shp")  # "../" means go up one level from the current folder or working directory
par(mar = c(2, 2, 2, 2))
plot(countries)
map.axes()  # add coordinates on the plot axes
# xmin, xmax, ymin, ymax coordinates of your region of interest:
my_window <- c(-22, 60, 25, 65)
plot(as(extent(my_window), "SpatialPolygons"), border = "red", lwd = 2, add = TRUE)  # check if it's where you want it on the map
# plot country borders within 'my_window' only:
plot(countries, xlim = my_window[1:2], ylim = my_window[3:4])
map.axes()


layers <- stack(list.files("./data/niche_model/climate/", pattern = "\\.tif$", full.names = TRUE))
plot(layers)

layers <- crop(layers, my_window)
plot(layers)
plot(layers[[1]], xlim = c(-1.631100, -0.998012), ylim = c(50.556800, 50.807402)) #check isle of wight
plot(layers[[1]], xlim = c(-2, 0), ylim = c(50, 51))

dat1 <- dat %>%
  mutate(latitude = as.character(latitude),
         longitude = as.character(longitude),
         latitude = if_else(population == "IOW2", "50.6", latitude),
         longitude = if_else(population == "IOW2", "-1.3", longitude)) %>% #have to change IOW pop. coordinates slightly as in raster with resolution I used they fall in the sea otherwise
  as.data.frame()
  
dat1 <- dat1 %>%  
  bind_cols(.,
            raster::extract(layers,
                            cbind(as.numeric(dat1$longitude), as.numeric(dat1$latitude)),
                            df = TRUE,
                            small = TRUE) %>%
              dplyr::select(-ID)) %>% 
  as.data.frame()

str(dat1)
View(dat1)

## Remove duplicates for each population (redundant)
### For populations with Nuc-Chl discordance:
### 1. I keep sample where there is match between Chl & Nuc if mism is in only one of two samples representing pop
### 2. I assign chloroplast group

pdn <- dat1 %>%
  filter(Plastome != Nuclear) %>% #View() #check which pops have mism, then only keep those that have mism for both samples
  filter(population == "14" | population == "Bro" | population == "Saf" | population == "Vil") %>%
  distinct(population, .keep_all = TRUE) %>%
  as.data.frame() #%>% View()

length(unique(dat$population)) #53: to check I do not drop any population by mistake

dat1 <- dat1 %>%
  group_by(population) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(seq_code != "L62") %>% #LIL because both nuclear and plastome assigned it to southwest
  filter(Plastome == Nuclear | n == 1) %>% 
  distinct(population, .keep_all = TRUE) %>% #View()
  dplyr::select(-n) %>%
  bind_rows(., pdn) %>%
  as.data.frame() #%>% View()

## Produce vector of dropped samples

rmvd.lb <- dat$seq_code[dat$seq_code %in% dat1$seq_code == FALSE]
rmvd <- c(rmvd.nlb, rmvd.lb, "L62")
rmvd


# ANALYSIS

## Create correlation matrix
set.seed(148)

mtx.clim <- cor(dat1[ , grep("^[WC_]", names(dat1), value=TRUE)]) #select 1 biologically meaningful variable per group when y < 2 (in paper it was 0.5, but they had many more datapoints)
heatmaply_cor(mtx.clim) #kept Tcv (bio4), Pcv (bio15), Tcoldestq (bio11), Pannual (bio12), Pdriestq (bio17)

## Check correlation & VIF between selected variables
mtx.clim2 <- cor(dat1[ , grep("bio4|bio15|bio11|bio12|bio17", names(dat1), value=TRUE)]) #bio15 and bio11 are still highly correlated (>0.8), but I keep them since they are biologically meaningful

vif(dat1[ , grep("bio4|bio15|bio11|bio17|bio12", names(dat1), value=TRUE)])
####WC_bio11_lonlat WC_bio12_lonlat WC_bio15_lonlat WC_bio17_lonlat  WC_bio4_lonlat 
####5.525559        2.690643        9.983649       11.698897        2.279497
####bio15 and bio17 vif are high (they are indeed similar vars, prec. seasonality and prec. driest quarter) > I remove bio17 with highest vif

vif(dat1[ , grep("bio4|bio15|bio11|bio12", names(dat1), value=TRUE)])
####WC_bio11_lonlat WC_bio12_lonlat WC_bio15_lonlat  WC_bio4_lonlat 
####5.291160        1.110127        3.946503        2.092559 
#### bio11 is at the edge (just above 5, used in paper as threshold to discard, for now I keep it)


## Run PCA with selected variables
set.seed(675)

nch.pca <- prcomp(dat1[ , grep("bio4|bio15|bio11|bio12", names(dat1), value=TRUE)], scale = TRUE)

screeplot(nch.pca, type = "l", npcs = 4, main = "Screeplot of the first 4 PCs")
abline(h = 1, col="red", lty = 5)
legend("topright", legend = c("Eigenvalue = 1"),
       col = c("red"), lty = 5, cex = 0.6)

cumpro <- cumsum(nch.pca$sdev^2 / sum(nch.pca$sdev^2))
plot(cumpro[0:4], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 4, col="blue", lty = 5)
abline(h = 0.9, col = "blue", lty = 5)
legend("topleft", legend = c("Cut-off @ PC4"),
       col = c("blue"), lty = 5, cex = 0.6)

pal.pc <- c("goldenrod", "darkblue", "skyblue", "firebrick")

fviz_pca_biplot(nch.pca,
                axes = c(1, 2),
                col.ind = dat1$Plastome, # Color by the quality of representation
                fill.ind = dat1$Plastome,
                geom = "point",
                pointsize = 3,
                stroke = 1.5,
                labelsize = 5,
                #addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                legend.title = "Plastome",
                palette = pal.pc,
                col.var = "gray42",
                label = c("var"),
                repel = TRUE     # Avoid text overlapping
                ) +
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  theme_bw(base_size = 24) +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1))

dat1 <- bind_cols(dat1, as.data.frame(nch.pca$x))

p12 <- ggplot(dat1, aes(x = PC1, y = PC2, col = Plastome, fill = Plastome, group = Plastome)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  bins = 5) +
  lims(x = c(-5, 5), y = c(-4, 4)) + 
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(nrow = 1), alpha = FALSE)

p13 <- ggplot(dat1, aes(x = PC1, y = PC3, col = Plastome, fill = Plastome, group = Plastome)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  bins = 5) +
  lims(x = c(-5, 5), y = c(-4, 4)) + 
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(nrow = 1), alpha = FALSE)

p14 <- ggplot(dat1, aes(x = PC1, y = PC4, col = Plastome, fill = Plastome, group = Plastome)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..),
                  bins = 5) +
  lims(x = c(-5, 5), y = c(-4, 4)) + 
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(nrow = 1), alpha = FALSE)



# Set up marginal histograms
x_hist1 <- ggplot(dat1, aes(x = PC1, fill = Plastome)) +
  geom_density(stat = "density", alpha = 0.4) +
  lims(x = c(-5, 5)) + 
  guides(fill = FALSE) +
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_void()

y_hist12 <- ggplot(dat1, aes(x = PC2, fill = Plastome)) +
  geom_density(stat = "density", alpha = 0.4) +
  lims(x = c(-4, 4)) + 
  guides(fill = FALSE) +
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_void() +
  coord_flip()

y_hist13 <- ggplot(dat1, aes(x = PC3, fill = Plastome)) +
  geom_density(stat = "density", alpha = 0.4) +
  lims(x = c(-4, 4)) +
  guides(fill = FALSE) +
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_void() +
  coord_flip()

y_hist14 <- ggplot(dat1, aes(x = PC4, fill = Plastome)) +
  geom_density(stat = "density", alpha = 0.4) +
  lims(x = c(-4, 4)) +
  guides(fill = FALSE) +
  scale_colour_manual(values = pal.pc) +
  scale_fill_manual(values = pal.pc) +
  theme_void() +
  coord_flip()

aligned_x_hist12 <- align_plots(x_hist1, p12, align = "v")[[1]]
aligned_y_hist12 <- align_plots(y_hist12, p12, align = "h")[[1]]

aligned_x_hist13 <- align_plots(x_hist1, p13, align = "v")[[1]]
aligned_y_hist13 <- align_plots(y_hist13, p13, align = "h")[[1]]

aligned_x_hist14 <- align_plots(x_hist1, p14, align = "v")[[1]]
aligned_y_hist14 <- align_plots(y_hist14, p14, align = "h")[[1]]


# Arrange plots
plot_grid(
  aligned_x_hist12
  , NULL
  , p12
  , aligned_y_hist12
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)

plot_grid(
  aligned_x_hist13
  , NULL
  , p13
  , aligned_y_hist13
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)

plot_grid(
  aligned_x_hist14
  , NULL
  , p14
  , aligned_y_hist14
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)


#prune beast phylogeny with ape (keep only one representative per group) OR prune it to keep one sample per population and assign PC values to all


## BEAST
tree <- read.nexus("./output/plastome/beast/all_plastomes_aligned_trimmed_supermx_beast_17NOempty.(time).MCMC.trees")
plot(tree, cex = 0.5)

right <- ladderize(tree, right = T)
plot(right)


### Ancestral state with all pops represented
right2 <- drop.tip(right, c(rmvd, "Linum_usitatissimum"))
plot(right2)

TipLabels(right2)
dat1$seq_code

TipLabels(right2) %in% dat1$seq_code

nch.mat1 <- setNames(dat1$PC1, dat1$seq_code)
nch.anc1 <- fastAnc(right2, nch.mat1)
anc.val1 <- contMap(right2, nch.mat1, plot = FALSE)
plot(anc.val1, type="phylogram", legend=0.7*max(nodeHeights(right2)))

nch.mat2 <- setNames(dat1$PC2, dat1$seq_code)
nch.anc2 <- fastAnc(right2, nch.mat2)
anc.val2 <- contMap(right2, nch.mat2, plot = FALSE)
plot(anc.val2, type="phylogram", legend=0.7*max(nodeHeights(right2)))

nch.mat3 <- setNames(dat1$PC3, dat1$seq_code)
nch.anc3 <- fastAnc(right2, nch.mat3)
anc.val3 <- contMap(right2, nch.mat3, plot = FALSE)
plot(anc.val3, type="phylogram", legend=0.7*max(nodeHeights(right2)))


### Ancestral state with one rep per group
tps <- TipLabels(right2)
kp <- c("L18", "L32", "L42", "L92")
rmvd2<- tps[!(tps %in% kp)]

right3 <- drop.tip(right2, rmvd2)
plot(right3)

dat2 <- dat1 %>%
  group_by(Plastome) %>%
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2),
            PC3 = mean(PC3),
            PC4 = mean(PC4)) %>%
  ungroup() %>%
  as.data.frame()

right3$tip.label <- c("CM", "South West", "North West", "East")

nch.mat12 <- setNames(dat2$PC1, dat2$Plastome)
nch.anc12 <- fastAnc(right3, nch.mat12)
anc.val12 <- contMap(right3, nch.mat12, plot = FALSE)
plot(anc.val12, type="phylogram", legend=0.7*max(nodeHeights(right3)))

nch.mat22 <- setNames(dat2$PC2, dat2$Plastome)
nch.anc22 <- fastAnc(right3, nch.mat22)
anc.val22 <- contMap(right3, nch.mat22, plot = FALSE)
plot(anc.val22, type="phylogram", legend=0.7*max(nodeHeights(right3)))

nch.mat32 <- setNames(dat2$PC3, dat2$Plastome)
nch.anc32 <- fastAnc(right3, nch.mat32)
anc.val32 <- contMap(right3, nch.mat32, plot = FALSE)
plot(anc.val32, type="phylogram", legend=0.7*max(nodeHeights(right3)))








