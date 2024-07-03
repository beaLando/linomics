library(tidyverse)
library(ggplot2)

# Gutaker sequence headers
gut.hds <- read.table("./data/Lb_TFL1/gutaker_headers.txt", sep = " ", header = FALSE)
str(gut.hds)

gut.hds <- gut.hds %>%
  rename(seq_code = "V1",
         accession = "V2") %>%
  mutate(seq_code = str_remove(seq_code, ">"),
         accession = str_remove(accession, "Lb")) %>%
  as.data.frame()

str(gut.hds)

# Gutaker infos
gut.info <- read.csv("./data/Lb_TFL1/gutaker_supplementary.csv")
str(gut.info)

gut.info <- gut.info %>%
  left_join(.,
            gut.hds,
            by = c("Donor_accession" = "accession")) %>%
  as.data.frame() #%>%
  #View()




# Portsmouth info
prt.info <- read.csv("./data/Lb_plastomes_raw/Consensus/populations_locations_1.csv")

toadd <- data.frame(
  spp = c("usitatissimum"),
  pop = c("cultivar"),
  seq_code = c("Linum_usitatissimum"),
  country = c(NA_character_),
  region = c("cultivated"),
  lat = c(NA),
  lon = c(NA),
  alt = c(NA),
  distance_coast = c(NA),
  location = c(NA_character_)
)

prt.info <- bind_rows(prt.info, toadd)

# All info

all.info <- gut.info %>%
  rename(spp = "Species",
         pop = "Donor_accession",
         seq_code = "seq_code",
         region = "region",
         country = "Country.of.origin",
         improvement_status = "Improvement_status",
         purpose = "Cultivation_purpose",
         lat = "Latitude",
         lon = "Longitiude",
         donor = "Donor",
         sample = "Sample") %>%
  select(-Collection_date, -Latitude_data, -Accession_name) %>%
  mutate(
    collection = rep("gutaker", nrow(.)),
    improvement_status = R.utils::decapitalize(improvement_status),
    purpose = R.utils::decapitalize(purpose),
    spp = if_else(spp == "bienne", "bienne", "usitatissimum")
    ) %>%
  mutate(purpose = if_else(is.na(purpose) == TRUE, improvement_status, purpose)) %>%
  bind_rows(.,
            prt.info %>%
              mutate(
                improvement_status = case_when(
                  pop == "cultivar" ~ "cultivar",
                  spp == "usitatissimum" ~ "unknown",
                  spp == "bienne" ~ "wild",
                  spp == "narbonense" ~ "wild"
                ),
                purpose = case_when(
                  improvement_status == "cultivar" ~ "oil",
                  spp == "usitatissimum" ~ "unknown",
                  spp == "bienne" ~ "wild",
                  spp == "narbonense" ~ "wild"
                ),
                donor = rep("Coll.", nrow(.)),
                sample = rep("collected", nrow(.)),
                collection = rep("portsmouth", nrow(.))
                )
              ) %>% #str()
  select(collection, sample, donor, spp, improvement_status, purpose, country, region, pop, seq_code, lat, lon, alt, LuTFL1, LuTFL2, RADseq) %>%
  as.data.frame()


# save file
write.csv(all.info, "./data/lblus_all_info.csv")






















