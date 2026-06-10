# cor.test(results$local_diversity, results$lat,method = 'spearman')

# Warning in cor.test.default(results$local_diversity, results$lat, method = "spearman") :
#   Cannot compute exact p-value with ties
# 
# Spearman's rank correlation rho
# data:  
# results$local_diversity
# results$lat
# S = 4875.3, p-value = 0.0001994
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6667807 

# Yann's NOTE:
# The radius around each (selfing, homozygous) individual was of 5 degrees (about 550km). I then extracted the average genetic distance between all pairs of individuals falling in this circle.
# The input files are attached, the piece of R code I used to generate the map is below:
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(geosphere)  # For distance calculations

# Load data (assumes columns: ID, lon, lat, and distance matrix)
#df <- read.table("Localities.txt",h=T)  # File with: ID, lon, lat
dist_matrix <- as.matrix(read.table("./data/Lb_genomes/yann_diversity_analysis/Distance_matrix_new.txt", row.names = 1))  # NxN matrix
colnames(dist_matrix)=rownames(dist_matrix)

df <- read.table("./data/Lb_genomes/yann_diversity_analysis/list_individuals_tokeep_selection.txt",h=T)  # File with: ID, lon, lat

df$ID <- gsub("L0([1-9])", "L\\1", df$seq_code)

# East and North-West lineage
df.enw=subset(df,df$Cluster_2>0.9 | df$Cluster_3>0.7)
dist_matrix.enw=dist_matrix[which(colnames(dist_matrix) %in% df.enw$ID),which(rownames(dist_matrix) %in% df.enw$ID)]

df.enw_ordered <- df.enw[match(rownames(dist_matrix.enw), df.enw$ID), ]
all(df.enw_ordered$ID == rownames(dist_matrix.enw))
df.enw=df.enw_ordered

# Define radius (in degrees of lat/lon)
radius <- 5# Adjust this to your preferred radius (e.g., 0.1 degrees)
map.polygon <- rworldmap::getMap(resolution = "low") ###Obtain maps from rworldmap

# Initialize results.enw
results.enw <- data.frame(lon = df.enw$lon, lat = df.enw$lat, local_diversity = NA)

# Compute local genetic diversity for each sample
for (i in 1:nrow(df.enw)) {
  focal_lon <- df.enw$lon[i]
  focal_lat <- df.enw$lat[i]

  # Compute distances between the focal point and all others
  dists <- distHaversine(c(focal_lon, focal_lat),cbind(df.enw$lon, df.enw$lat)) / 1000  # Convert meters to km

  # Find individuals within the given radius (converted to km)
  nearby <- which(dists <= radius * 111)  # 1 degree ≈ 111 km

  # Compute mean genetic distance within local area
  if (length(nearby) > 1) {  # Ensure we have multiple individuals (at least 2)
    local_distances <- dist_matrix.enw[i, nearby]
    results.enw$local_diversity[i] <- mean(local_distances[local_distances > 0], na.rm = TRUE)
  }
}

# South-West Lineage
df.sw=subset(df,df$Cluster_1>0.9)
dist_matrix.sw=dist_matrix[which(colnames(dist_matrix) %in% df.sw$ID),which(rownames(dist_matrix) %in% df.sw$ID)]

df.sw_ordered <- df.sw[match(rownames(dist_matrix.sw), df.sw$ID), ]
all(df.sw_ordered$ID == rownames(dist_matrix.sw))
df.sw=df.sw_ordered

# Initialize results.sw
results.sw <- data.frame(lon = df.sw$lon, lat = df.sw$lat, local_diversity = NA)

# Compute local genetic diversity for each sample
for (i in 1:nrow(df.sw)) {
  focal_lon <- df.sw$lon[i]
  focal_lat <- df.sw$lat[i]
  
  # Compute distances between the focal point and all others
  dists <- distHaversine(c(focal_lon, focal_lat),cbind(df.sw$lon, df.sw$lat)) / 1000  # Convert meters to km
  
  # Find individuals within the given radius (converted to km)
  nearby <- which(dists <= radius * 111)  # 1 degree ≈ 111 km
  
  # Compute mean genetic distance within local area
  if (length(nearby) > 1) {  # Ensure we have multiple individuals (at least 2)
    local_distances <- dist_matrix.sw[i, nearby]
    results.sw$local_diversity[i] <- mean(local_distances[local_distances > 0], na.rm = TRUE)
  }
}

results <- results.enw %>%
  bind_cols(.,
            df.enw) %>%
  mutate(group = "NW + East Lineage") %>%
  filter(!is.na(group)) %>%
  droplevels() %>%
  bind_rows(.,
            results.sw %>%
            bind_cols(.,
            df.sw) %>%
              mutate(group = "SW Lineage") %>%
              filter(!is.na(group)) %>%
              droplevels()) %>%
  as.data.frame() #lat and lon redundant across dfs, but I leave it because I can check with View() if there is correspondence

results.enw <- droplevels(subset(results, group == "NW + East Lineage"))
results.sw <- droplevels(subset(results, group == "SW Lineage"))

# Plot the map
ggplot() +
  # Add map outline (optional, requires "map.polygon" with coastline data)
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group), color = "black", alpha = 0.5) +
  
  # Plot local diversity as colored circles
  geom_point(data = results.enw, aes(x = lon...1, y = lat...2, fill = local_diversity),
             shape = 21, size = 6, color = "black") +
  scale_fill_gradient2(low = "lightgray", high = "blue", midpoint=0.06, na.value = "gray") +
  labs(fill = "Diversity - East and NW Lineage") +

  # start a new scale
  new_scale_fill() +
  
  # apply the black-grey gradient to group 2
  geom_point(data = results.sw, aes(x = lon...1, y = lat...2, fill = local_diversity), 
             shape = 21, size = 6, color = "black") +
  scale_fill_gradient2(low = "lightgray", high = "red", midpoint=0.09, na.value = "gray") +
  labs(fill = "Diversity - SW Lineage") +
  
  # Plot actual sample locations as black dots
  geom_point(data = results, aes(x = lon...1, y = lat...2), color = "black", size = 2) +

  # Labels & Theme
  xlim(-20,40) +
  ylim(25,57) +
  coord_equal() +
  xlab("Longitude") + ylab("Latitude") + 
  theme_bw(base_size = 18) +
  theme(legend.position = "top")
