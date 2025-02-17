library(tidyverse)
library(magrittr)
library(terra)
library(tidyterra)
library(giscoR)
library(sf)

nuts0 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "0")

germany <- gisco_get_nuts(
  year = "2021",
  country = 'DE',
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "0")
nuts0_cropped <-  st_crop(nuts0,
                          xmin = -11, 
                          ymin =  34,
                          ymax = 72 ,
                          xmax =  34 )


corine_all <- rast('./data/raster_data/corine/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif')
pastures <- corine_all == 'Pastures'
binary <- as.int(pastures)
#pasture_10kmsq <- aggregate(binary, fact=10, fun = 'mean')

#binary_polgons <- as.polygons(binary)
#binary_sf 
#pasture_10kmsq_df <- as.data.frame(binary) %>% as_tibble()
#pasture_10kmsq_sf <- binary %>% 
#  as.polygons() %>%
 # st_as_sf() %>%
 # st_transform( crs = 4326)

test <- st_intersection(pasture_10kmsq_sf, nuts0)


sextent <- terra::ext(germany)
binary_cropped <- terra::crop(binary, sextent)
crs(binary) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
binary_masked <- terra::mask(binary, vect(germany))

plot(binary_masked)
ggplot() + 
  geom_sf(data = nuts0_cropped) + 
 # coord_sf()
  #geom_sf(data = pasture_10kmsq_sf, aes(fill  = LABEL3), colour = NA) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) 



geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) + <- 

st_as_sf(pasture_10kmsq_df, coords = c('x', 'y'), crs = st_crs(corine_all))


y <- classify(pastures, cbind(TRUE, 1))

pastures[, , T] <- 1
pastures[, , F] <- 0

pastures_df <- as.data.frame(pastures, xy = TRUE, na.rm  = TRUE)
pastures_true <- pastures_df %>% filter(layer == TRUE)

true_sf <- st_as_sf(pastures_true, coords = c('x', 'y'), crs = st_crs(corine_all))
geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) + <- 