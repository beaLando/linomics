# LIBRARIES
library(tidyverse)
library(ggplot2)
library(leaflet)
library(GGally)
library(parzer)
library(xlsx)

# FIXING COORDINATES

## Import pop info with wrong coordinates for Gutaker
pop.info <- xlsx::read.xlsx("./data/populations_locations_1.xlsx", sheetName = "pop_info_plastome")
str(pop.info)

sppCol <- colorFactor(palette = 'RdYlGn', pop.info$spp)

leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addCircleMarkers(
    data = pop.info,
    lng = ~as.numeric(lon),
    lat = ~as.numeric(lat),
    popup = ~as.character(pop), #~as.character(Pop)
    color = ~sppCol(spp)) 

# NB: After checking map, some of Gutaker's L. bienne samples seem off (fall in sea)
## It could be that I made a mistake in DMS coordinate translation to decimal degrees.
## Check Gutaker's thesis and retrieve coordinates for those samples again (also do L. usitatissimum while you are at it).

### Gutaker's samples current coords
# usitatissimum	M072	L70	Russia	55	103	NA
# usitatissimum	M073	L71	Kazakhstan	48.005284	66.9045434	NA
# usitatissimum	M074	L72	Ukraine	48.383022	31.1828699	NA
# usitatissimum	M075	L73	Ukraine	48.383022	31.1828699	NA
# usitatissimum	M076	L74	Kazakhstan	48.005284	66.9045434	NA
# bienne	W65	L75	Croatia	43.06805556	16.11888889	358
# bienne	W70	L76	Greece	39.06694444	20.03583333	281
# bienne	W72	L77	Greece	38.15083333	21.15083333	1
# bienne	W74	L78	Greece	39.05166667	21.035	946
# bienne	W76	L79	Greece	39.11861111	21.13583333	621
# bienne	W77	L80	Greece	39.13416667	21	835
# bienne	W81	L81	Greece	40.01722222	23.1	2
# bienne	W82	L82	Greece	40.08388889	23.11861111	129
# bienne	W85	L83	Bulgaria	41.08444444	23.05111111	591
# bienne	W88	L84	Bulgaria	42.74444444	27.8375	89
# bienne	W94	L85	Croatia	45.13333333	18.06722222	149
# bienne	W95	L86	Croatia	45.03583333	16.10194444	268
# bienne	W96	L87	Croatia	45.08361111	15.03416667	180

### Gutaker's samples DMS coords from thesis
new.coords <- data.frame(
  pop = c("M072", "M073", "M074", "M075", "M076", #Lus from PGR Canada - 30852, 30841, 30844, 30846, 30842
          "W65", "W70", "W72", "W74", 
          "W76", "W77", "W81", "W82", 
          "W85", "W88", "W94", "W95", #NB I've also noticed that, while Allaby shared W88 with us as "L. bienne" in Gutaker's thesis this sample is marked as "?" as species > Indeed, this sample was not aligning well with all other L. bienne samples in phylogeny and we had to exclude it.
          "W96"),
  seq_code = paste0("L", 70:87),
  latlon.dms = c(NA, NA, NA, NA, NA,
                 "N43°54'555 E016°27'880", "N39°54'514 E020°22'296", "N38°59'438 E021°09'233", "N39°03'265 E021°52'565",
                 "N39°47'476 E021°38'599", "N39°58'435 E021°30'800", "N40°21'625 E023°56'005", "N40°35'028 E023°47'074",
                 "N41°35'840 E023°43'843", "N42°43'100 E027°45'315", "N45°08'308 E018°14'326", "N45°22'092 E016°16'175",
                 "N45°25'711 E015°22'231"),
  alt = c(NA, NA, NA, NA, NA,
          358, 281, 1, 946, 
          621, 835, 2, 129,
          591, 89, 149, 268,
          180),
  country = c("Rostov, Russia", "Kazakhstan", "Ukraine", "Ukraine", "Kazakhstan",
              "Croatia", "Greece", "Greece", "Greece", 
              "Greece", "Greece", "Greece", "Greece",
              "Bulgaria", "Bulgaria", "Croatia", "Croatia",
              "Croatia"))


## Retrieve decimal coordinates with parzer
new.coords <- new.coords %>%
  separate(latlon.dms, into = c("lat.dms", "lon.dms"), sep = " ") %>%
  mutate(lat = parse_lat(lat.dms),
         lon = parse_lon(lon.dms)) %>%
  filter(grepl("^W", pop)) %>% # I keep only bienne samples because usitatissimum looks ok
  mutate(spp = if_else(pop == "W88", "unknown", "bienne")) %>%
  as.data.frame() #%>% View()


## Fix info dataframe & replot for checks
pop.info <- pop.info %>%
  filter(!(pop %in% c(new.coords$pop))) %>%
  bind_rows(.,
            new.coords %>%
              mutate_all(as.character) %>%
              select(-c(lat.dms, lon.dms))) %>%
  as.data.frame()

sppCol <- colorFactor(palette = 'RdYlGn', pop.info$spp)

leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addCircleMarkers(
    data = pop.info,
    lng = ~as.numeric(lon),
    lat = ~as.numeric(lat),
    popup = ~as.character(pop), #~as.character(Pop)
    color = ~sppCol(spp)) 

View(pop.info)


## Overwrite excel file (then manually correct other sheet with coordinates in same file)
wb <- loadWorkbook("./data/populations_locations_1.xlsx")
removeSheet(wb, sheetName = "pop_info_plastome")
yourSheet <- createSheet(wb, sheetName="pop_info_plastome")
addDataFrame(pop.info, yourSheet)
saveWorkbook(wb, "./data/populations_locations_1_coordFixed.xlsx")
