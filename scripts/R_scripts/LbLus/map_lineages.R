# Re-coding seq_code
# 2022-2-22

library(readxl)
library(tidyverse)

# Read table 2
T2 <- read_excel("G:/Shared drives/PHDthesis/5_Phylogeography/v2/Figures/Table2_Lbienne_Lusitatissimum_samples.xlsx") 

glimpse(T2)

# Plastome 
T2 <- T2 %>% 
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

# Nuclear
T2 <- T2 %>% 
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


# TFL1:
# 88 & 84 : Out groups
# 74-65 : Other
# 44-69 : North Med
# 73-76 : Cultivated
# 87-29 : West
# 61-22 : South West 1
# 30-40 : South West 2

T2 <- T2 %>% 
  mutate(.,
         TFL1 = recode(seq_code,
                       "L88" = "Outgroups",
                       "L74" = "Other",
                       "L55" = "Other",
                       "L94" = "Other",
                       "L10" = "Other",
                       "L93" = "Other",
                       "L49" = "Other",
                       "L47" = "Other",
                       "L54" = "Other",
                       "L11" = "Other",
                       "L60" = "Other",
                       "L53" = "Other",
                       "L80" = "Other",
                       "L92" = "Other",
                       "L81" = "Other",
                       "L79" = "Other",
                       "L96" = "Other",
                       "L82" = "Other",
                       "L12" = "Other",
                       "L13" = "Other",
                       "L95" = "Other",
                       "L91" = "Other",
                       "L64" = "Other",
                       "L65" = "Other",
                       "L84" = "Outgroups",
                       "L44" = "North Med",
                       "L56" = "North Med",
                       "L59" = "North Med",
                       "L03" = "North Med",
                       "L01" = "North Med",
                       "L02" = "North Med",
                       "L18" = "North Med",
                       "L68" = "North Med",
                       "L41" = "North Med",
                       "L69" = "North Med",
                       "L73" = "Cultivated",
                       "L70" = "Cultivated",
                       "L72" = "Cultivated",
                       "L71" = "Cultivated",
                       "L78" = "Cultivated",
                       "L43" = "Cultivated",
                       "L77" = "Cultivated",
                       "L83" = "Cultivated",
                       "L76" = "Cultivated",
                       "L87" = "West",
                       "L63" = "West",
                       "L58" = "West",
                       "L05" = "West",
                       "L67" = "West",
                       "L04" = "West",
                       "L48" = "West",
                       "L06" = "West",
                       "L42" = "West",
                       "L51" = "West",
                       "L66" = "West",
                       "L45" = "West",
                       "L75" = "West",
                       "L86" = "West",
                       "L52" = "West",
                       "L07" = "West",
                       "L50" = "West",
                       "L26" = "West",
                       "L37" = "West",
                       "L15" = "West",
                       "L46" = "West",
                       "L14" = "West",
                       "L85" = "West",
                       "L89" = "West",
                       "L27" = "West",
                       "L62" = "West",
                       "L90" = "West",
                       "L25" = "West",
                       "L17" = "West",
                       "L57" = "West",
                       "L19" = "West",
                       "L32" = "West",
                       "L34" = "West",
                       "L35" = "West",
                       "L36" = "West",
                       "L33" = "West",
                       "L20" = "West",
                       "L21" = "West",
                       "L28" = "West",
                       "L29" = "West",
                       "L61" = "South West 1",
                       "L08" = "South West 1",
                       "L09" = "South West 1",
                       "L23" = "South West 1",
                       "L22" = "South West 1",
                       "L30" = "South West 2",
                       "L31" = "South West 2",
                       "L39" = "South West 2",
                       "L38" = "South West 2",
                       "L40" = "South West 2",
                       "Lus" = NA_character_))



# Plotting ----------------------------------------------------------------

library(ggplot2)
library(ggnewscale)

table(T2$Plastome)
table(T2$Nuclear)
table(T2$TFL1)

table(T2$TFL1)

T2 <- T2 %>% 
  filter(Plastome != "Outgroups") %>% 
  filter(seq_code != "Lus") %>% 
  mutate(.,
         TFL1_nr = recode(TFL1,
                          "West" = 1,
                          "South West 1" = 1,
                          "North Med" = 2,
                          "Cultivated" = 3,
                          "South West 2" = 4,
                          "Other" = 5)) %>% 
  mutate_at(.,
            vars(longitude, latitude),
            as.numeric) %>% 
  filter(spp != "usitatissimum")


# PAPER

T3 <- T2 %>%
  bind_cols(.,
            raster::extract(maxent_proj,
                            cbind(as.numeric(T2$longitude), as.numeric(T2$latitude)),
                            df = TRUE) %>%
              rename(pres.suit = "layer") %>%
              select(-ID)) %>%
  bind_cols(.,
            raster::extract(maxent_proj_hol, 
                            cbind(as.numeric(T2$longitude), as.numeric(T2$latitude)),
                            df = TRUE) %>%
              rename(holo.suit = "layer") %>%
              select(-ID)) %>%
  bind_cols(.,
            raster::extract(maxent_proj_lgm, 
                            cbind(as.numeric(T2$longitude), as.numeric(T2$latitude)),
                            df = TRUE) %>%
              rename(lgm.suit = "layer") %>%
              select(-ID)) %>%
  bind_cols(.,
            raster::extract(maxent_proj_lig,
                            cbind(as.numeric(T2$longitude), as.numeric(T2$latitude)),
                            df = TRUE) %>%
              rename(lig.suit = "layer") %>%
              select(-ID)) %>%
  bind_cols(.,
            raster::extract(maxent_proj_plio, 
                            cbind(as.numeric(T2$longitude), as.numeric(T2$latitude)),
                            df = TRUE) %>%
              rename(plio.suit = "layer") %>%
              select(-ID)) %>%
  as.data.frame()


T3 %>%
  filter(Plastome != "Outgroups") %>%
  filter(spp != "usitatissimum") %>%
  pivot_longer(pres.suit:plio.suit,
               names_to = "epoch",
               values_to = "suitability") %>%
  mutate(epoch = gsub(".suit", "", epoch)) %>% #View()
  mutate(epoch = factor(epoch, levels = c("plio", "lig", "lgm", "holo", "pres")),
         Plastome = factor(Plastome, levels = c("CM", "South West", "North West", "East"))) %>% #View()#str()
  group_by(epoch, Plastome, population) %>%
  summarise(suitability = mean(suitability, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = epoch, y = suitability, col = Plastome)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1),
             alpha = 0.3) +
  # geom_line(aes(group = seq_code),
  #           position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1),
  #           alpha = 0.1,
  #           size = 0.1) +
  stat_summary(fun = "mean", 
               aes(group = Plastome),
               geom = "line", 
               position = position_dodge(width = 0.5),
               size = 0.5) +
  stat_summary(fun.data = "mean_se", 
               position = position_dodge(width = 0.5),
               size = 1) +
  scale_colour_manual(values = c("gold", "red2", "deepskyblue1", "dodgerblue4")) +
  #facet_wrap(~Plastome, nrow = 2) +
  theme_minimal(base_size = 15)







# POSTER

maxent_proj.df <- as.data.frame(as(maxent_proj, "SpatialPixelsDataFrame"))
maxent_proj_lgm.df <- as.data.frame(as(maxent_proj_lgm, "SpatialPixelsDataFrame"))

mapl <- ggplot() +  
  geom_tile(data = maxent_proj_lgm.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  scale_fill_gradient(low = "black", high = "white") +
  coord_equal() +
  theme_minimal() +
  theme(legend.position="right")

mapp <- ggplot() +  
  geom_tile(data = maxent_proj.df, aes(x = x, y = y, fill = layer), alpha = 0.8) + 
  scale_fill_gradient(low = "black", high = "white") +
  coord_equal() +
  theme_minimal() +
  theme(legend.position="bottom")

m.app <- mapp + 
  #coord_quickmap(xlim = c(-20, 36), ylim = c(29, 54)) +
  new_scale_fill() +
  geom_point(data = T2, size = 2, stroke = 2, shape = 21, 
             aes(x = longitude, y = latitude,
                 #shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  #scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "dodgerblue4", 'deepskyblue1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "dodgerblue4", 'deepskyblue1', "red2")) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("right"))
m.app



#THESIS

world_map <- map_data("world")

mapa <- ggplot(legend = FALSE) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray") +
  theme_minimal()

# Generic map
png(filename = "Generic_map.png", units = "in", res = 400, width = 4, height = 5)
m.all <- mapa + 
  # coord_quickmap(xlim = c(-20, 35), ylim = c(29, 55)) +
  coord_quickmap(xlim = c(-20, 36), ylim = c(29, 54)) +
  geom_point(data = T2, size = 1.25, 
             aes(x = longitude, y = latitude)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c("none"))
m.all
dev.off()

# Generic map with only plastome and nuclear (no tfl1)
png(filename = "Generic_map_cols.png", units = "in", res = 400, width = 5, height = 4)
m.all <- mapa + 
  coord_quickmap(xlim = c(-20, 36), ylim = c(29, 54)) +
  geom_point(data = T2, size = 2, stroke = 2, shape = 21, 
             aes(x = longitude, y = latitude,
                 #shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  #scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "dodgerblue4", 'deepskyblue1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "dodgerblue4", 'deepskyblue1', "red2")) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("none"))
m.all
dev.off()  


# Insets
png(filename = "Inset_A_map.png", units = "in", res = 400, width = 5, height = 4)
m.all <- mapa + 
  coord_quickmap(xlim = c(-18, 9), ylim = c(28, 43.8)) +
  # coord_quickmap(xlim = c(-14, 35), ylim = c(28, 62)) +
  geom_point(data = T2, size = 2, stroke = 2, 
             aes(x = longitude, y = latitude,
                 shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c("none"))
m.all
dev.off()  


png(filename = "Inset_B_map.png", units = "in", res = 400, width = 6, height = 5)
m.all <- mapa + 
  coord_quickmap(xlim = c(9, 36), ylim = c(33, 46)) +
  # coord_quickmap(xlim = c(-14, 35), ylim = c(28, 62)) +
  geom_point(data = T2, size = 2, stroke = 2, 
             aes(x = longitude, y = latitude,
                 shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("none"))

m.all
dev.off()  


png(filename = "Inset_C_map.png", units = "in", res = 400, width = 6, height = 5)
m.all <- mapa + 
  coord_quickmap(xlim = c(-5, 6), ylim = c(44.5, 54)) +
  # coord_quickmap(xlim = c(-14, 35), ylim = c(28, 62)) +
  geom_point(data = T2, size = 2, stroke = 2, 
             aes(x = longitude, y = latitude,
                 shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c("none"))

m.all
dev.off()  


png(filename = "Inset_D_map.png", units = "in", res = 400, width = 3, height = 3)
m.all <- mapa + 
  coord_quickmap(xlim = c(-7, -2.3), ylim = c(35, 39)) +
  geom_point(data = T2, size = 2, stroke = 2, 
             aes(x = longitude, y = latitude,
                 shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c("none"))

m.all
dev.off()  

png(filename = "Inset_E_map.png", units = "in", res = 400, width = 3, height = 3)
m.all <- mapa + 
  coord_quickmap(xlim = c(-2.5, 0), ylim = c(52, 54)) +
  geom_point(data = T2, size = 2, stroke = 2, 
             aes(x = longitude, y = latitude,
                 shape = as.factor(TFL1_nr),
                 color = Nuclear,
                 fill = Plastome,
             )) +
  scale_shape_manual(name = "TFL1", values = c(21, 22, 23, 24, 25))+
  scale_color_manual(name = "Lineages", values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  scale_fill_manual(name = "Lineages",values = c("gold", "steelblue", 'paleturquoise1', "red2")) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c("none"))

m.all
dev.off()  


