# LIBS
library(tidyverse)
library(sf)
library(raster)
library(fuzzySim)
library(maps)

# IMPORT DATA
sppd <- read.csv("./data/niche_model/species_occurrence/linumspp_1700today_geoclean.csv")
layers <- raster::stack(list.files(path = "./data/niche_model/climate/wc2.1_10m_bio/", pattern = "*_crpd.tif$", full.names = TRUE))
countries <- st_read("./data/niche_model/countries/world_countries.shp")  # "../" means go up one level from the current folder or working directory

# SELECT Lb and Lus

spp <- c("linum bienne")
sppd <- sppd %>% filter(species %in% spp) %>% dplyr::select(-c(X.1, X)) %>% as.data.frame()

sppd.lb <- sppd %>% filter(species %in% spp[1]) %>% as.data.frame()

remove(sppd)

# DELIMIT THE MODELLING REGION ####

# convert species occurrences to a spatial object:
occurrences_spatial.lb <- sppd.lb

names(occurrences_spatial.lb)

coordinates(occurrences_spatial.lb) <- sppd.lb[ , c("decimalLongitude", "decimalLatitude")]
raster::crs(occurrences_spatial.lb) <- "+proj=longlat"

plot(occurrences_spatial.lb)
plot(countries, border = "blue", col = "transparent", lwd = 2, add = TRUE)
map.axes()


# select the modelling region, e.g. using the countries where this species has occurrence points (which means the species was surveyed in those countries):
countries_with_points.lb <- subset(countries, ADMIN == "France" | ADMIN == "Spain" | ADMIN == "Portugal" | ADMIN == "Andorra")

plot(occurrences_spatial.lb)
plot(countries_with_points.lb, border = "red", col = "transparent", lwd = 2, add = TRUE)
map.axes()

# judge if some countries are visibly insufficiently surveyed (in this dataset) for this species; compare with occurrence data from other sources, e.g. https://www.iucnredlist.org, national atlases, data papers
# also, only the mainland should be included in the modelling region, as on islands species may be limited by dispersal
# so, assuming we can't complement the dataset with occurrences from other sources (which we could!), let's select mainland Spain for modelling:

# create a unique polygon identifier for 'countries_with_points' and add it as labels to the map, to see which polygon(s) we want to select:
countries_with_points.lb$my_id <- 1:nrow(countries_with_points.lb)
raster::text(countries_with_points.lb, labels = countries_with_points.lb$my_id, col = "blue", font = 2, halo = TRUE)
# select only the desired polygon(s) for the modelling region (polygon 1 FOR THE EXAMPLE DATA - CHANGE AS APPROPRIATE!!!):
mod_region.lb <- subset(countries_with_points.lb, my_id %in% c(14, 5, 32, 24, 23)) #c(3, 9, 11)
plot(mod_region.lb, border = "green", lwd = 4, add = TRUE)  # check that it is the desired polygon(s)

# add species points and filter them with the modelling region:
plot(occurrences_spatial.lb, col = "grey", add = TRUE)
occurrences_spatial.lb <- occurrences_spatial.lb[sf::as_Spatial(mod_region.lb),]
plot(occurrences_spatial.lb, col = "darkblue", add = TRUE)

# IF YOU USED A LIMITED WINDOW OF COORDINATES to download the occurrence data, you need to intersect or crop with that too:
my_window <- c(-24, 79, 22, 70)

mod_region.lb <- crop(sf::as_Spatial(mod_region.lb), extent(my_window))
plot(mod_region.lb, border = "darkgreen", lwd = 3, add = TRUE)

# now import and cut (crop + mask) the variable maps to the extent of the modelling region:
##LB
layers_mod.lb <- mask(crop(layers, mod_region.lb), mod_region.lb)
names(layers_mod.lb)
plot(layers_mod.lb[[1]])
plot(countries, border = "darkgray", col = "transparent", lwd = 2, add = TRUE)
plot(mod_region.lb, add = TRUE)
plot(occurrences_spatial.lb, add = TRUE)


# SET THE APPROPRIATE SPATIAL RESOLUTION ####
# closely inspect your species data vs. the size of the variables' pixels:

plot(layers_mod.lb[[1]], xlim = c(0.5, 2), ylim = c(44, 46))
points(occurrences_spatial.lb, cex = 0.5)
plot(layers_mod.lb[[1]], xlim = c(-8, -6), ylim = c(37, 39))
points(occurrences_spatial.lb, cex = 0.5)

# # plot within different x/y limits if necessary to see if presence point resolution matches pixel resolution (i.e., if you don't have evenly spaced points with pixels in between -- see lecture "04_predictor_variables.pdf")
# # notice that, in the example data, pixels have approximately the same spatial resolution as the presence points, but they don't match exactly (i.e., points are not at the centroids of these pixels); ideally, you should find the original grid over which (most of) these presences were sampled, and extract the raster values to that grid
# 
# # if necessary (which is not the case for the example data), you can aggregate the layers, to e.g. a 5-times coarser resolution (choose the 'fact' value that best matches your presence data resolution to your variables' resolution):
# layers_aggr.lb <- raster::aggregate(layers_mod.lb, fact = 2, fun = mean)
# 
# res(layers_aggr.lb)
# 
# plot(layers_aggr.lb[[1]], xlim = range(occurrences_spatial.lb$decimalLongitude), ylim = range(occurrences_spatial.lb$decimalLatitude))
# points(occurrences_spatial.lb, pch = ".")

# save the modelling layers to a folder on disk:
saveRDS(layers_mod.lb, paste0("./data/niche_model/climate_mod_", spp[1], ".rds"))

# now make a dataframe of the species occurrence data gridded to the resolution of the raster variables
# i.e., one row per pixel with the values of the variables and the presence/absence of species records:

head(sppd.lb)

## LB
gridded_data.lb <- gridRecords(rst = layers_mod.lb, pres.coords = sppd.lb[ , c("decimalLongitude", "decimalLatitude")])
head(gridded_data.lb)

nrow(gridded_data.lb)  # should be the same number as:
sum(!is.na(getValues(layers_mod.lb[[1]])))

names(gridded_data.lb)

# plot the gridded records:
plot(layers_mod.lb[[1]])
# plot the absences (pixels without presence records):
points(gridded_data.lb[gridded_data.lb[ , "presence"] == 0, c("x", "y")], col = "red", cex = 0.5)
# plot the presences (pixels with presence records):
points(gridded_data.lb[gridded_data.lb[ , "presence"] == 1, c("x", "y")], col = "blue", cex = 0.7)

# plot within a narrower coordinate range to see closer:
plot(occurrences_spatial.lb)
map.axes()
plot(layers_mod.lb[[1]], xlim = c(-8, -4), ylim = c(36, 39))
points(gridded_data.lb[gridded_data.lb[ , "presence"] == 0, c("x", "y")], col = "red", cex = 0.5)
points(gridded_data.lb[gridded_data.lb[ , "presence"] == 1, c("x", "y")], col = "blue", pch = 20)


# save the modelling dataframe to a .csv file on disk:
write.csv(gridded_data.lb, paste0("./data/niche_model/dat_", spp[1], ".csv"), row.names = FALSE)
