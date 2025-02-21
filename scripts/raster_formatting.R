####################################################################################################
####################################################################################################
## Script name: Format Rasters
##
## Purpose of script: Extract and format environmental and ecological raster layers for analysis in
## Seraphim (Dellicour et al. 2016)
##
## Date created: 2025-02-20
##
##
########################################## SYSTEM OPTIONS ##########################################


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(terra)
library(tidyterra)
library(giscoR)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(spatstat)


# User functions


############################################## DATA ################################################
europe_sf <- bind_rows(
  
  # GISCO resgions
  gisco_get_nuts(
    year = "2021",
    epsg = "4326", #WGS84 projection
    resolution = "03", #1:10million
    nuts_level = "0"),  
  
  # Non GISCO regions
  ne_countries(scale = 10, 
               country = c('bosnia and herzegovina', 'kosovo', 'andorra'),
               returnclass = "sf") %>%
    dplyr::select(name_en,
                  iso_a2_eh,
                  geometry) %>%
    rename(CNTR_CODE = iso_a2_eh,
           NAME_LATN = name_en) %>%
    mutate(NUTS_ID = CNTR_CODE)) %>%
  
  # Filter
  st_crop(ont, xmin=-11, xmax=34, ymin=34, ymax=72)


# Untransformed CORINE data
#corine_all <- rast('./data/raster_data/corine/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif')

# Transformed CORINE data
corine_all_masked <- rast('./data/raster_data/corine/u2018_clc2018_v2020_20u1_raster100m/separate_layers/masked_cropped.tiff')

############################################## MAIN ################################################
### Mask and refine CORINE Data ##

europe <- vect(europe_sf)

# Change CORINE projection and coodinate system
# corine_all_latlong <- project(corine_all, "+proj=longlat epsg:4326")
# remove(corine_all)

# Only keep CORINE data within Europe mask
# corine_all_cropped <- crop(corine_all_latlong, ext(-11, 34, 34, 72))
# remove(corine_all_latlong)
# corine_all_masked <- mask(corine_all_cropped, europe)
# remove(corine_all_cropped)

# Extract single layer
pastures <- corine_all_masked == 'Pastures'
#remove(corine_all_latlong)

corine_classes <- levels(corine_all_masked)[[1]]
classes_to_include <-  corine_classes %>% pull(LABEL3) %>% magrittr::extract(c(21:26, 28:43))


SeparateAndSave <- function(raster, layer){
  separate_layer <- raster == layer
  
  filename <- str_to_lower(layer) %>%
    str_replace_all(., ' ' , '_' ) %>%
    paste0('./data/raster_data/corine/u2018_clc2018_v2020_20u1_raster100m/separate_layers/', ., '.tiff')
  
  writeRaster(separate_layer, filename)
  
}

SaveAsPointsDF <- function(rasterfile){
  r <- rast(rasterfile)
  
  
  
  
}
sep_layers <- lapply(classes_to_include, SeparateAndSave, raster = corine_all_masked)  
### To plot CORINE data using geom_spatraster

# Represent as presence/absence
binary <- as.int(pastures)
#remove(pastures)

ggplot() +
  geom_spatraster(data = urban,  maxcell = 1e+06) +
  scale_fill_viridis_c(na.value = "grey90") + 
  geom_sf(colour = 'white', data = europe_sf, linewidth = 0.2, alpha = 0.01) + 
  coord_sf(expand = FALSE)


### Extract CORINE data as points and plot as smoothed density 
pasture_points <- as.data.frame(pastures, xy=TRUE,  na.rm=TRUE)
urban_points <- as.data.frame(urban, xy=TRUE,  na.rm=TRUE)

ggplot() +
  #geom_spatraster(data = urban,  maxcell = 1e+06) +
  geom_sf(data = urban_points_sf, size = 0.5, color = "black")+
  #scale_fill_viridis_c(na.value = "grey90") + 
  geom_sf(colour = 'white', data = europe_sf, linewidth = 0.2, alpha = 0.01) + 
  coord_sf(expand = FALSE)


urban_points_ppp <- as.ppp(urban_points_sf_2$geometry)


europe_sf_2 <- st_transform(europe_sf, st_crs("ESRI:102118"))
ga_crs <- st_crs("ESRI:102118")
urban_points_sf_2 <- st_transform(urban_points_sf, st_crs("ESRI:102118"))

# Convert the church coordinates to a ppp object with a built-in window
urban_points_ppp <- as.ppp(urban_points_sf_2$geometry, W = as.owin(europe_sf_2))

# Create a stars object (whatever that is) of the density of church locations
density_churches_stars <- stars::st_as_stars(density(ga_churches_ppp, dimyx = 300))

# Convert the stars object to an sf object so it's normal and plottable again
ga_churches_density <- st_as_sf(density_churches_stars) %>%
  st_set_crs(ga_crs)
############################################## WRITE ###############################################


############################################## END #################################################
####################################################################################################
####################################################################################################













# make coastline mask
coast_vector <- europe |> 
  # transfrom to crs of raster
  st_transform(crs = terra::crs(binary)) |> 
  # convert to spatVector
  vect()  


binary_points = as.points(binary)

#pasture_10kmsq <- aggregate(binary, fact=10, fun = 'mean')

#binary_polgons <- as.polygons(binary)
#binary_sf 
#pasture_10kmsq_df <- as.data.frame(binary) %>% as_tibble()
#pasture_10kmsq_sf <- binary %>% 
#  as.polygons() %>%
# st_as_sf() %>%
# st_transform( crs = 4326)




test <- st_intersection(pasture_10kmsq_sf, nuts0)


ggplot(d) + 
  geom_sf()

sextent <- terra::ext(germany)
binary_cropped <- terra::crop(binary, sextent)
crs(binary) <- "EPSG:4326"
binary_masked <- terra::mask(binary, vect(germany))

plot(binary_masked)
ggplot() + 
  #geom_sf(data = nuts0) + 
  geom_spatraster(data = binary) + 
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