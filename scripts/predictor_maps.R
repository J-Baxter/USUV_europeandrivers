# Geospatial predictors
library(sf)
library(rnaturalearth)
library(tidyverse)

# Import data
test <- readRDS('./spatial_data/Predictors_tep_20231014.Rdata')


# import europe data
#map <- ne_countries(returnclass = "sf", scale = 'medium')

country <- read_sf('../spatial_data/gadm_410-levels.gpkg', layer = 'ADM_0')
subdivision_one <- read_sf('../spatial_data/gadm_410-levels.gpkg', layer = 'ADM_1')
subdivision_two <- read_sf('../spatial_data/gadm_410-levels.gpkg', layer = 'ADM_2')


ggplot(data = country) +
  geom_sf(aes()) +
  geom_point(data = test, mapping = aes(x = lon, y = lat, color = tmp))+
  theme_bw()