################################################################################
## Script Name:        Format Rasters
## Purpose:            Extract and format environmental and ecological raster 
##                     layers for analysis in Seraphim (Dellicour et al. 2016). 
##                     For consistency, we grid all rasters using the ETRS89 10km
##                     grid.
## Author:             James Baxter
## Date Created:       2025-09-19
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(terra)
library(tidyterra)
library(sf)
library(stars)


FilterCorine <- function(raster, selected_layers){
  out <- raster %in% selected_layers
  return(out)
}


AggregateCorine <- function(raster, grid_vector){
  # Turn logical raster into presence/absence binary
  raster <- as.int(raster)
  out <- exact_extract(raster,
                       grid_vector,
                       fun = 'mean') %>% 
    bind_cols(grid_vector, 'value' = .) #zonal(raster, ETRS89_10_v, fun = mean, na.rm = T)
  
  return(out)
  
}


SaveLayer <- function(sf_dataframe, raster_name){
  filename <- str_to_lower(raster_name) %>%
    str_replace_all(., ' ' , '_' ) %>%
    str_replace_all(., '[:punct:]' , '' ) %>%
    paste0('./2025Jun24/raster_data/formatted/', .)
  
  # as raster
  raster <- st_rasterize(sf_dataframe %>% dplyr::select(value, geometry))
  write_stars(raster, paste0(filename, '.tif'))
  
  # as shapefile
  st_write(sf_dataframe, paste0(filename, '.shp'))
  
  cat(paste0('Saved ', raster_name, ' successfully.'))
}


################################### DATA #######################################
# Read and inspect data

# Import ETRS89 grid
ETRS89_10 <- read_sf('./2025Jun24/raster_data/ebba/ebba2_etrs10x10_v1/ebba2_etrs10x10_v1.shp') 
#ETRS89_10_v <- vect(ETRS89_10)

# Untransformed CORINE data
corine_all <- rast('./2025Jun24/raster_data/corine/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif')

# CORINE classes
corine_legend <- read_csv('./2025Jun24/raster_data/clc_legend.csv')

# EBBA2 data
ebba_csv <- read_delim('./2025Jun24/raster_data/ebba/ebba2_data_model_10km.csv', delim = ';')

################################### MAIN #######################################
# Main analysis or transformation steps

# 1. CORINE
# Only keep CORINE data within Europe mask
corine_all_masked <- mask(corine_all, ETRS89_10)
remove(corine_all)

# Agregate corine classes by level two
level_two_dict <- corine_legend %>%
  drop_na() %>%
  split(~ LABEL2) %>%
  lapply(., function(x) x %>% pull(LABEL3)) %>%
  as.list()

#corine_grouped <- list()
#corine_grouped[[1]] <- FilterCorine(corine_all_masked, level_two_dict[[1]])
corine_grouped <- lapply(level_two_dict, FilterCorine, raster = corine_all_masked) %>%
  setNames(names(level_two_dict))

#corine_agregated[[1]] <- AggregateCorine(corine_grouped[[1]], ETRS89_10)
corine_agregated <- lapply(corine_grouped, AggregateCorine, grid_vector = ETRS89_10) %>%
  setNames(names(level_two_dict))

mapply(SaveLayer, 
       corine_agregated, 
       names(level_two_dict))


# 2. EBBA
birds <- ebba_csv %>%
  pull(birdlife_scientific_name) %>% 
  unique()

ebba_gridded <- ebba_csv %>%
  rename(value = 'Probability of occurrence') %>%
  split(~ birdlife_scientific_name) %>%
  map(~ ETRS89_10 %>% left_join(.x)) 

mapply(SaveLayer, 
       ebba_gridded, 
       names(ebba_gridded))

# 3. FAO 

# 4. Elevation

# 5. Human Population

# 6. ERA5 (Monthly)

# 7. Bioclimatic
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################



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
l


# User functions


############################################## DATA ################################################


############################################## MAIN ################################################
### Mask and refine CORINE Data ##

europe <- vect(europe_sf)



# Extract single layer
pastures <- rast('./2025Jun24/raster_data/pastures.tiff')

#remove(corine_all_latlong)



SaveAsPointsDF <- function(rasterfile){
  r <- rast(rasterfile)
  
  
  
  
}
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

pasture_10kmsq <- aggregate(binary, fact=10, fun = 'mean')

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