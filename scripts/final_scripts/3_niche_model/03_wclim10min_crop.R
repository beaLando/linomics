path <- "./data/niche_model/climate/wc2.1_10m_bio/" #NB: 30s grid is 0.08333333 resolution // 10m grid is 0.166667 resolution
fls <- list.files(path)

stck <- stack(as.list(paste0(path, fls)))

#area to extract for (lon, lon, lat, lat)
ext <- c(-25, 80, 21, 71)

stck.cropped <- crop(stck, ext)  
plot(stck.cropped[[1]])

names(stck.cropped) <- gsub("wc2.1_10m_", "", names(stck.cropped))

#write data to files
writeRaster(stck.cropped, filename = paste0(path, paste0(names(stck.cropped), "_crpd")), bylayer = TRUE, format = "GTiff")


