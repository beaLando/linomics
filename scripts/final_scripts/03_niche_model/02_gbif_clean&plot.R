# Libs
library(tidyverse)

library(rgbif)
library(countrycode)
library(CoordinateCleaner)

library(maps)

library(ggplot2)
library(ggfortify)
library(patchwork)


## Load data & get columns useful for cleaning occurrences

linum.gbif0 <- read.csv("./data/niche_model/species_occurrence/linumspp_1700today.csv")

linum.gbif <- linum.gbif0 %>%
  mutate_all(as.character) %>% 
  mutate_all(~iconv(.x,"WINDOWS-1252","UTF-8")) %>% # some variables contain unreadable characters
  mutate_all(tolower) %>%
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude),
         coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>%
  as.data.frame()

# Get count of obs for Linum species in EuropeXAngiosperm353
spp <- c("bienne")
spp <- paste0("linum ", spp)

linum.gbif0.summary <- linum.gbif %>%
  filter(species %in% spp) %>% 
  droplevels() %>% #str()
group_by(species) %>%
  summarise(n.tot = n()) %>%
  ungroup() %>%
  as.data.frame() #%>% str()


## Quick map the occurrence data

length(unique(linum.gbif$species))

map("world", 
    xlim = range(linum.gbif$decimalLongitude), 
    ylim = range(linum.gbif$decimalLatitude))  # if the map doesn't appear right at first, run this command again

points(linum.gbif[ , c("decimalLongitude", "decimalLatitude")], col = "darkgray", pch = ".")
points(linum.gbif[linum.gbif$species == "linum bienne" , c("decimalLongitude", "decimalLatitude")], col = "red", pch = ".")


## Quick clean (also check Marcia's code) #https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/

### Get records which refer to any Linum spp. being absent from location
unique(linum.gbif$occurrenceStatus)
absnt.nms <- c("absent", "ausente", "assente", "absente", "abwesend", "afwezig", "frånvarande") #could add more

length(which(linum.gbif$individualCount == 0 | linum.gbif$occurrenceStatus %in% absnt.nms)) #indCount refers to how many Linum spp. individuals were recorded during survey at location // alternative is just present/absent under occStatus
#there is 1 absence

linum.gbif <- linum.gbif %>%
  mutate(to.rm = if_else(
    individualCount == 0 | occurrenceStatus %in% absnt.nms, "remove_abs", "keep")
    ) %>%
  as.data.frame()


### Get records with impossible/unlikely/... coords

nrow(linum.gbif) #24255

linum.sub <- linum.gbif %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude)) %>%
  mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
         countryCode = countrycode(countryCode,
                                   origin =  'iso2c',
                                   destination = 'iso3c')) %>%
  cc_val() %>%
  cc_equ() %>%
  cc_cap() %>%
  cc_cen() %>%
  #cc_coun(iso3 = "countryCode") %>%
  cc_sea() %>%
  cc_zero() %>%
  cc_outl() %>%
  cc_dupl() %>%
  as.data.frame()

nrow(linum.sub) #15626

linum.flagged <- linum.gbif %>%
  filter(!(X %in% linum.sub$X)) %>%
  as.data.frame()

nrow(linum.flagged) #8629


### Plot by to.rm

map("world", 
    xlim = range(linum.gbif$decimalLongitude), 
    ylim = range(linum.gbif$decimalLatitude))  # if the map doesn't appear right at first, run this command again

points(linum.sub[ , c("decimalLongitude", "decimalLatitude")], col = "blue", pch = ".")
points(linum.flagged[ , c("decimalLongitude", "decimalLatitude")], col = "red", pch = ".")


## Save data with clean coords

linum.sub %>%
  select(c(decimalLongitude, decimalLatitude, gbifID, species, scientificName, eventDate)) %>% 
  arrange(species, eventDate) %>% 
  View()

write.csv(linum.sub, "./data/niche_model/species_occurrence/linumspp_1700today_geoclean.csv") #could be improved













