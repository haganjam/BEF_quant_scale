
# map the study sites
library(sp)
library(sf)
library(raster)
library(dplyr)
library(readr)
library(ggplot2)

# load the Swedish coastline
swe_coast <- st_read("data/GIS_layers/Sweden_coast/land_skagerrak_kattegat.shp")
projection(swe_coast)

# convert to wgs84 coordinate lat lon reference system
swe_coast <- st_transform(swe_coast, crs = "+proj=longlat +datum=WGS84 +no_defs")

# plot the Swedish coastline
# plot(swe_coast)

# get a map of just the Tjarno region
tja_coast <- st_crop(swe_coast, xmin = 11, xmax = 11.3, ymin = 58, ymax = 59)
plot(tja_coast, col = "white", bg = "white", main = NULL)

# get the coordinates of the sampling points
cs2 <- read_csv("data/case_study_2/ResearchBox 843/Data/site_data.csv")
head(cs2)

# get the latitude-longitude coordinates
cs2 <- dplyr::select(cs2, cluster_id, dec_lat, dec_lon)
head(cs2)

# get the unique cluster ids
cs2 <- 
  cs2 %>%
  group_by(cluster_id) %>%
  summarise(dec_lat = mean(dec_lat),
            dec_lon = mean(dec_lon))

# remove the cluster F because we did not use it
cs2 <- 
  cs2 %>%
  filter(cluster_id != "F")

# rename the columns
cs2_pts <- 
  cs2 %>%
  dplyr::select(dec_lon, dec_lat)
names(cs2_pts) <- c("longitude", "latitude")

# Convert data frame to sf object
cs2_pts <- st_as_sf(x = cs2_pts, coords = c("longitude", "latitude"), crs = st_crs(swe_coast))

# plot the points
plot(tja_coast, col = "white")
points(cs2_pts, col = "red")

ggplot() +
  geom_sf(data = tja_coast) +
  geom_sf(data = cs2_pts)




