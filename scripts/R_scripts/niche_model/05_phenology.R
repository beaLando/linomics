occurrences_spatial@data$year_add <- rep("2020", nrow(occurrences_spatial@data))

occurrences_spatial@data$date_obs <- paste0(as.character(occurrences_spatial@data$day), 
                                            "-", 
                                            as.character(occurrences_spatial@data$month), 
                                            "-", 
                                            occurrences_spatial@data$year_add)

occurrences_spatial@data$date_obs <- as.Date(occurrences_spatial@data$date_obs, format = "%d-%m-%Y")
occurrences_spatial@data$date_baseline <- rep(as.Date("2020-01-01"), nrow(occurrences_spatial@data))

occurrences_spatial@data$phenology_obs <- as.integer((occurrences_spatial@data$date_obs - occurrences_spatial@data$date_baseline))

ex.df <- as.data.frame(raster::extract(gam_proj_fav, occurrences_spatial, cellnumbers=T))
ex.df <- cbind(ex.df, xyFromCell(gam_proj_fav, ex.df[ ,1]))

names(ex.df) <- c("gam_cellID", "gam_fav", "lon_cell", "lat_cell")

occurrences_spatial@data <- bind_cols(occurrences_spatial@data, ex.df)

pheno_obs <- occurrences_spatial@data %>%
  group_by(gam_cellID) %>%
  summarise(
    lon_cell = unique(lon_cell),
    lat_cell = unique(lat_cell),
    gam_fav = unique(gam_fav),
    mean_pheno = as.integer(mean(phenology_obs, na.rm = TRUE)),
    sd_pheno = sd(phenology_obs, na.rm = TRUE),
    mean_yr = as.integer(mean(as.numeric(year))),
    sd_yr = sd(as.numeric(year), na.rm = TRUE),
    n = n()
  ) %>%
  filter(is.na(lon_cell) == FALSE) %>%
  filter(is.na(lat_cell) == FALSE) %>%
  filter(is.na(gam_fav) == FALSE) %>%
  filter(mean_pheno > 0) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(gam_fav_q = cut(gam_fav,
                         breaks = quantile(gam_fav, c(0, 0.25, 0.5, 0.75, 1)),
                         labels = c("1-25", "26-50", "51-75", "76-100"),
                         na.rm = TRUE)
         ) %>%
  filter(is.na(gam_fav_q) == FALSE) %>%
  droplevels() %>%
  as.data.frame()

coordinates(pheno_obs) <- pheno_obs[ , c("lon_cell", "lat_cell")]
crs(pheno_obs) <- "+proj=longlat"
par(mfrow = c(1, 1))
plot(pheno_obs)

ggplot(pheno_obs@data, aes(x = lat_cell, y = gam_fav, col = mean_pheno)) +
  geom_point() +
  #geom_smooth() +
  theme_bw()

ggplot(pheno_obs@data, aes(x = gam_fav_q, y = mean_pheno)) +
  geom_jitter(aes(col = lat_cell)) +
  geom_boxplot(alpha = 0.01) +
  #geom_smooth() +
  theme_bw()


plot(gam_proj_fav, col = clrs, breaks = brks, main = "GAM")
points(occurrences_spatial[ , c("decimalLongitude", "decimalLatitude")], cex = 0.01, pch = 20, col = "black")  # compare e.g.with the range map of this species at https://www.iucnredlist.org to assess if the distribution is well represented
points(pheno_obs[ , c("lon_cell", "lat_cell")], cex = 0.01, pch = 20, col = "red")

































