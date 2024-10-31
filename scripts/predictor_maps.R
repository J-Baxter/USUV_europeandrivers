# Geospatial predictors
library(sf)
library(rnaturalearth)
library(tidyverse)

# Import data
test <- readRDS('./spatial_data/Predictors_tep_20231014.Rdata')


# import europe data
#map <- ne_countries(returnclass = "sf", scale = 'medium')

nuts0 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "0")

nuts1 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "1")

nuts2 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "2")

nuts3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "3")

all_europe_sf <- bind_rows(nuts0, nongisco_shapefile)

gdp <- test %>%
  dplyr::select(lat,lon,GDP)

gdp_sf <- st_as_sf(gdp, coords = c('lon', 'lat'), crs = 4326)

europe_map <- ne_countries(continent = 'Europe', scale= 'medium', returnclass = 'sf')
grid_size <- 10000

europe_bbox <- st_transform(europe_map, crs = st_crs(3857))
europe_grid <- st_make_grid(europe_bbox, cellsize = grid_size, what = 'polygons')
europe_grid <- st_sf(geometry = europe_grid)

europe_grid <- st_intersection(europe_grid, st_transform(europe_map, 3857))

gdp_sf <- st_transform(gdp_sf, 3857)

grid_gdp <- st_join(europe_grid, gdp_sf, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(gdp = sum(GDP, na.rm = T))

ggplot(grid_gdp)+
  geom_sf(aes(fill = gdp))+
  theme_minimal()

ggplot(data = all_europe_sf) +
  geom_sf( colour = 'black', fill = 'lightgrey') + 
  geom_point(data = test, mapping = aes(x = lon, y = lat, color = popcount))+
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  theme_void() + 
  theme(legend.position = 'none')