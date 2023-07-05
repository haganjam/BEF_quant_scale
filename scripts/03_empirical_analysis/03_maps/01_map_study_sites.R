
# map the study sites
library(sp)
library(sf)
library(raster)
library(dplyr)
library(readr)
library(ggplot2)
library(ggspatial)

# load plotting theme
source("scripts/Function_plotting_theme.R")

# case study 1

# load the British coastline
br_coast <- st_read("data/GIS_layers/GBR_adm/GBR_adm0.shp")
projection(br_coast)

# convert to wgs84 coordinate lat lon reference system
br_coast <- st_transform(br_coast, crs = "+proj=longlat +datum=WGS84 +no_defs")
showDefault(br_coast)

# get the Plymouth region
ply_coast <- st_crop(br_coast, xmin = -6, xmax = -3, ymin = 49.8, ymax = 50.7)

# check the crop
ggplot() +
  geom_sf(data = ply_coast)

# get the coordinates of Challaborough and Kingsand
cs1_pts <- 
  tibble(site = c("CB", "KS"),
         dec_lon = c(-3.9, -4.195),
         dec_lat = c(50.2833, 50.333333))

# convert the coordinates to an sf object
cs1_pts <- st_as_sf(x = cs1_pts, coords = c("dec_lon", "dec_lat"), crs = st_crs(br_coast))

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 9, type = "continuous")

# plot the map
p1 <- 
  ggplot() +
  geom_sf(data = ply_coast) +
  geom_sf(data = cs1_pts, mapping = aes(fill = cluster_id),
          size = 2.5, fill = col_pal[5], shape = 21, colour = "black", stroke = 0.5) +
  annotation_scale(location = "tr", width_hint = 0.2, height = unit(0.175, "cm")) +
  coord_sf(xlim = c(-4.3, -3.8), ylim = c(50.20, 50.43), expand = FALSE) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(0.75, "cm"),
                         width = unit(0.75, "cm"),
                         pad_x = unit(0.05, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering(text_size = 9, line_width = 0.1) ) +
  theme_void()
plot(p1)

ggsave("figures/MAIN_fig_2a.svg", p1, units = "cm",
       width = 6.5, height = 4)


# case study 2

# load the Swedish coastline
swe_coast <- st_read("data/GIS_layers/Sweden_coast/land_skagerrak_kattegat.shp")
projection(swe_coast)

# convert to wgs84 coordinate lat lon reference system
swe_coast <- st_transform(swe_coast, crs = "+proj=longlat +datum=WGS84 +no_defs")

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
  summarise(dec_lat = first(dec_lat),
            dec_lon = first(dec_lon))

# remove the cluster F because we did not use it
cs2 <- 
  cs2 %>%
  filter(cluster_id != "F")

# rename the columns
cs2_pts <- 
  cs2 %>%
  dplyr::select(cluster_id, dec_lon, dec_lat)

# Convert data frame to sf object
cs2_pts <- st_as_sf(x = cs2_pts, coords = c("dec_lon", "dec_lat"), crs = st_crs(swe_coast))

p2 <- 
  ggplot() +
  geom_sf(data = tja_coast) +
  geom_sf(data = cs2_pts, mapping = aes(fill = cluster_id),
          size = 2.75, shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = col_pal) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(xlim = c(11.09, 11.20), ylim = c(58.845, 58.91), expand = FALSE) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         pad_x = unit(0.05, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering(text_size = 9, line_width = 0.1) ) +
  guides(fill = guide_legend(override.aes = list(size = 2),
                             nrow = 1,
                             label.position = "bottom",
                             label.vjust=7.5)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.001, 'cm'),
        legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),
                                   size = 7))
plot(p2)

ggsave("figures/MAIN_fig_2b.svg", p2, units = "cm",
       width = 7, height = 8.5)

### END
