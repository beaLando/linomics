library(tidyverse)
library(rgbif)
library(ISOcodes)
library(lingtypology)

scientific_name <- c("Linum bienne", "Linum usitatissimum")

vernacular_names1 <- bind_rows(name_lookup(query = scientific_name[1])$names) %>%
  mutate(scientificName = rep(scientific_name[1], nrow(.))) %>%
  as.data.frame()

vernacular_names2 <- bind_rows(name_lookup(query = scientific_name[2])$names) %>%
  mutate(scientificName = rep(scientific_name[2], nrow(.))) %>%
  as.data.frame()

vernacular_names <- bind_rows(vernacular_names1,
                              vernacular_names2) %>%
  rename(language.iso = language) %>%
  left_join(.,
            ISO_639_2 %>%
              select(Alpha_3_T, Name),
            by = c("language.iso" = "Alpha_3_T")) %>%
  rename(language = Name) %>%
  select(scientificName, language, language.iso, vernacularName) %>%
  #mutate(country = country.lang(language),
         #lon = long.lang(language.iso),
         #lat = lat.lang(language.iso)) %>%
  #select(scientificName, country, language, language.iso, lon, lat, vernacularName) %>%
  as.data.frame()

str(vernacular_names)
unique(vernacular_names$language)
