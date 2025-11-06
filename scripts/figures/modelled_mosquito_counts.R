################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-11-06
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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(giscoR)


################################### DATA #######################################
# Read and inspect data
metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')

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

nuts_all <- bind_rows(nuts0,
                      nuts1,
                      nuts2)

vector_net <- read_sf('./spatial_data/VectornetMAPforMOODjan21.shp', 
                      crs = st_crs(nuts0)) %>%
  st_make_valid()

#vector_net <- st_transform(vector_net, st_crs(vector_net))

culex_abundance_files <- list.files('./culex_models/statistical_models/Abundace_model_predictions',
                                    full.names = T)

################################### MAIN #######################################
# Main analysis or transformation steps
# Main analysis or transformation steps
# Pre-format vectornet polygons
vect_id <- as_tibble(vector_net[-1083,]) %>%
  rowid_to_column() %>%
  st_as_sf()

# Prepare predictions from abundance data
culex_abundance <- lapply(culex_abundance_files, read_csv) %>%
  setNames(gsub('\\.\\/culex_models\\/statistical_models\\/Abundace_model_predictions\\/culex_|\\.csv', 
                '',
                culex_abundance_files)) %>%
  bind_rows(., .id = 'mm_yyyy') %>%
  separate_wider_delim(., mm_yyyy, '_', names = c('month', 'year')) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) 

# Simple average (mean) over all month-year combinations, aggregated by 
# polygons
culex_abundance_by_vectornet <- aggregate(culex_abundance %>% 
                                            dplyr::select(prediction, geometry),
                                          vect_id, 
                                          mean) %>%
  rowid_to_column() %>%
  mutate(prediction = log10(prediction))


# In some cases (e.g London), there is insufficient resolution to aggregate NUTS3
# polygons. In these cases, we estimate the polygon value as the median of it's
# non-NA neighbours
interpolate_missing_polygons <- culex_abundance_by_vectornet %>% 
  filter(is.na(prediction)) %>%
  mutate(eurostat_polygon = map(geometry, ~ {
    which(st_intersects(., culex_abundance_by_vectornet, sparse = FALSE)[1,])
  })) %>%
  st_drop_geometry() %>%
  
  # Expand dataframe so that there is 1 row per eurostate polygon
  unnest(eurostat_polygon) %>%
  filter(rowid != eurostat_polygon) %>%
  
  # join eurostat polygons, with aggregated culex abundance. Selecting best
  # resolution by default
  left_join(culex_abundance_by_vectornet,
            by = join_by('eurostat_polygon' == 'rowid')) %>%
  st_drop_geometry() %>% 
  summarise(prediction = median(prediction.y, na.rm = T), .by = rowid) %>%
  drop_na()

culex_abundance_by_vectornet %<>%
  as_tibble() %>% 
  rows_update(., interpolate_missing_polygons, by = 'rowid') %>%
  st_as_sf() 


################################### OUTPUT #####################################
# Save output files, plots, or results
ggplot(culex_abundance_by_vectornet,
       aes(fill = prediction)
       
) + 
  geom_sf(linewidth = 0) + 
  geom_sf(data = nuts1, colour = 'grey75', fill = '#FF000000') + 
  geom_sf(data = nuts0, colour = 'black', fill = '#FF000000') + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdYlGn',
                       'Predicted Mosquito Count per Month',
                       #transform = 'log10', 
                       direction = -1, , 
                       na.value="lightgrey") +
  ##facet_grid(rows = vars(year), cols = vars(month)) + 
  theme_void() +
  theme(legend.position = 'bottom')

ggsave('./2025Jun24/plots/modelled_mosquito_counts.jpeg',
       height = 15,
       width = 15,
       units = 'cm',
       dpi =360)

#################################### END #######################################
################################################################################
