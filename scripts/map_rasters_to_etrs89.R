################################################################################
## Script Name:        Format Rasters
## Purpose:            Extract and format environmental and ecological raster 
##                     layers for analysis in Seraphim (Dellicour et al. 2016). 
##                     For consistency, we grid all rasters to the 10 km x 10 km 
##                     ETRS89 grid, using the ETRS89-extended / LAEA Europe 
##                     (EPSG:3035) projection.
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
library(exactextractr)
library(giscoR)

FilterCorine <- function(raster, selected_layers){
  out <- raster %in% selected_layers
  return(out)
}


AggregateLayer <- function(raster, grid_vector, transform = TRUE){
  if(isTRUE(transform)){
    # Turn logical raster into presence/absence binary
    raster <- as.int(raster)
  }
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
  sf::write_sf(sf_dataframe, 
               paste0(filename, '.shp'))
  
  cat(paste0('Saved ', raster_name, ' successfully.\n'))
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

# FAO data
fao_list <- list.files('./2025Jun24/raster_data/fao',
                       pattern = 'tif',
                       full.names = T)
fao_rast <- lapply(fao_list, rast)

# EU-DEM data
elevation_original <- rast('./2025Jun24/raster_data/EU_DEM_mosaic_1000K/eudem_dem_3035_europe.tif')

#ERA5 monthly
era5_all <- rast('./2025Jun24/raster_data/era5_data.grib')


################################### MAIN #######################################
# Main analysis or transformation steps
# Set bounding box
bounding_box <- st_bbox(c(xmin = 2500000, xmax = 6700000, 
                          ymin = 1200000, ymax = 5500000), 
                        crs = st_crs(ETRS89_10)) %>%
  st_as_sfc()

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

corine_grouped <- lapply(level_two_dict, 
                         FilterCorine, 
                         raster = corine_all_masked) %>%
  setNames(names(level_two_dict))

corine_agregated <- lapply(corine_grouped,
                           AggregateLayer, 
                           grid_vector = ETRS89_10) %>%
  lapply(., function(x) x[st_intersects(x, bounding_box, sparse = FALSE),]) %>%
  setNames(names(level_two_dict))


# 2. EBBA
birds <- ebba_csv %>%
  pull(birdlife_scientific_name) %>% 
  unique()

ebba_gridded <- ebba_csv %>%
  rename(value = 'Probability of occurrence') %>%
  split(~ birdlife_scientific_name) %>%
  map(~ ETRS89_10 %>% left_join(.x)) %>%
  lapply(., function(x) x[st_intersects(x, bounding_box, sparse = FALSE),])


# 3. FAO 
# Rescale and mask to europe
fao_reprojected <- fao_rast %>%
  lapply(., function(x) {
    x <- project(x, 'EPSG:3035')
    return(x)}) 

fao_gridded <- lapply(fao_reprojected,
                      AggregateLayer, 
                      grid_vector = ETRS89_10)%>%
  lapply(., function(x) x[st_intersects(x, bounding_box, sparse = FALSE),])

fao_names <- paste0('GLW4-2020_', c('chicken', 'cattle', 'goat', 'swine', 'sheep'))


# 4. Elevation
elevation_gridded <- elevation_original %>%
  AggregateLayer(grid_vector = ETRS89_10) %>%
  .[st_intersects(., bounding_box, sparse = FALSE),] 


# 5. Eurostat Human Population
log10_population <- gisco_get_grid(
  resolution = "10",
  spatialtype = c("REGION")) %>%
  .[st_intersects(., bounding_box, sparse = FALSE),] %>%
  dplyr::select(TOT_P_2018, geom) %>%
  rename(geometry = geom) %>%
  mutate(value = log10(TOT_P_2018)) %>%
  dplyr::select(-TOT_P_2018)


# 6. ERA5 (Monthly)
era5_names <- names(era5_all) %>%
  str_extract(., "(?<=;)[^\\(\\[\\{]+") %>%
  str_trim() %>%
  gsub(".*\\/", '', .) %>%
  paste(time(era5_all) %>% format("%Y%b"),
        .,
        sep = '_')

names(era5_all) <- era5_names

era5_list <- lapply(names(era5_all), function(nm) era5_all[[nm]])

era5_gridded <- lapply(era5_list, 
                       AggregateLayer,
                       grid_vector = ETRS89_10, 
                       transform = FALSE) %>%
  lapply(., function(x) x[st_intersects(x, bounding_box, sparse = FALSE),]) %>%
  setNames(era5_names)

# 7. Bioclimatic


################################### OUTPUT #####################################
# Save output files, plots, or results
# 1. CORINE
mapply(SaveLayer, 
       corine_agregated, 
       names(level_two_dict))

# 2. EBBA
mapply(SaveLayer, 
       ebba_gridded, 
       names(ebba_gridded))

# 3. FAO
mapply(SaveLayer, 
       fao_gridded, 
       fao_names)

# 4. EU-DEM
SaveLayer(elevation_gridded, 'elevation')

# 5. Eurostat Human Population
SaveLayer(log10_population, 'log10_humanpopulation')

# 6. ERA5 - monthly
mapply(SaveLayer, 
       era5_gridded, 
       era5_names)

#fao_gridded[[2]] %>%
#ggplot() + 
#geom_sf(aes(fill = value), colour = '#ff000000') + 
#coord_sf(datum = sf::st_crs(ETRS89_10),
#ylim = c(1200000, 5500000),
#xlim = c(2500000, 6700000), 
#expand = FALSE) 
#scale_fill_continuous(trans = 'log10')
#################################### END #######################################
################################################################################