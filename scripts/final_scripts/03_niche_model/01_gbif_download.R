# Libs
library(tidyverse)

library(rgbif)
library(CoordinateCleaner)

library(maps)

library(ggplot2)
library(ggfortify)
library(patchwork)


## Import gbif occurrences (can validate phenology with Pl@ntNet pics + Herbaria)

spps <- c("Linum bienne")
coords <- data.frame(Lat = c("22", "70"), Lon = c("-24", "79"))
years <- c("1700, 1899", #divided into smaller chunks with passing of time because otherwise download is too big
           "1900, 1950", "1951, 1970", "1971, 1990", "1991, 1995", "1996, 2000", 
           "2001, 2005", "2006, 2010", "2011, 2015", "2016, 2020", "2021, 2023") #, 2024-2025: left this out because maybe data not all there yet

linum.gbif <- list()

for(i in 1:length(years)){
  
  dt.i<- occ_data(
    scientificName = spps, 
    hasCoordinate = TRUE, 
    decimalLongitude = paste0(coords$Lon, collapse = ", "), decimalLatitude = paste0(coords$Lat, collapse = ", "),
    year = years[i],
    limit = 500000)
  
  linum.gbif[[i]] <- as.data.frame(dt.i$data)
  
} #FOR CITATION: GBIF.org (03 July 2025) GBIF Occurrence Download  https://doi.org/10.15468/dl.nhbuhr


## Clean Linum df from doubles based on records IDs, location and date of record, species names
#linum.gbif.previous <- read.csv("./data/linumspp_1700today.csv")

linum.gbif <- bind_rows(linum.gbif) %>% #str()
  select(-networkKeys) %>%
  mutate_all(as.character) %>% #str()
  # bind_rows(., #if dates are updated and more records are downloaded
  #           linum.gbif.previous %>%
  #             mutate_all(as.character)) %>%
  distinct(key, datasetKey, gbifID, occurrenceID,
           scientificName, 
           decimalLatitude, decimalLongitude, eventDate, #eventDate is when occurrence was observed (https://www.gbif.org/data-quality-requirements-occurrences#dcEventDate), while dateIdentified should be the date on which species was IDed
           .keep_all = TRUE)


## Save Linum df
write.csv(as.data.frame(linum.gbif), "./data/niche_model/species_occurrence/linumspp_1700today.csv")






