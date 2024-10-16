# paper figures
library(terra)
library(sf)
library(rnaturalearth)
library(giscoR)
library(eurostat)
library(tidyverse)


# Initialise map of europe subdivisions (NUTS level 2)
# Data obtained from Eurostat GISCO (the Geographic Information System of the COmmission) via API 
# package giscoR

nuts2 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "10", #1:10million
  nuts_level = "2"
)

ggplot(nuts2) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE)



# Import data
lu_spatialdata <- readRDS('./spatial_data/Predictors_tep_20231014.Rdata')

#test raster
test_df <- lu_spatialdata %>% 
  dplyr::select(c(lon, lat, Deciduous_broad))


ggplot(nuts2) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE) + 
  geom_raster(aes(y = lat, x = lon, fill = Deciduous_broad), data = test_df)


coordinates(test_df) <- ~ lon + lat
r <- raster(extent(test_df), resolution = 0.1)
crs(r) <- CRS("+proj=longlat +datum=WGS84")

raster_test <- rasterize(test_df, r, field = 'Deciduous_broad', fun = mean)

library(raster)

raster::extract()


ggplot(data = country) +
  geom_sf(aes()) +
  geom_point(data = test, mapping = aes(x = lon, y = lat, color = tmp))+
  theme_bw()