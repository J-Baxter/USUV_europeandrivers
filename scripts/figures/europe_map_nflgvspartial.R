################################################################################
## Script Name:        NUTS-3 Map for sequence distribution
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-07-25
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


GetUncertaintyWeightings <- function(data, geo_data, nflg_only = FALSE){
  out <- data %>% 
    filter(is_europe == '1') %>%
    {
      if (nflg_only) {
        filter(., generegion_nflg == '1') 
      } else {
        .
      }
    } %>%
    filter(location_precision %in% c('nuts0', 'nuts1', 'nuts2')) %>%
    
    # format unique IDs
    mutate(seq_id = coalesce(sequence_accession, sequence_isolate) %>%
             gsub('\\/', '_', .)) %>%
    dplyr::select(seq_id, ends_with('_id')) %>%
    pivot_longer(cols = -1, names_to = 'nuts_level', values_to = 'nuts_id') %>%
    
    # filter out unnecessary levels
    drop_na() %>%
    mutate(nuts_level = gsub('nuts|_id', '', nuts_level) %>% as.integer(.)) %>%
    slice_max(nuts_level, by = seq_id) %>%
    #filter(seq_id =='Blackbird_London_2022_a') %>%
    
    # Join current polygon
    left_join(geo_data %>% dplyr::select(NUTS_ID, geometry),
              by = join_by('nuts_id' == 'NUTS_ID')) %>%
    
    # list all vectornet polygons within current polygon
    mutate(eurostat_polygon = map(geometry, ~ {
      which(st_overlaps(vect_id$geometry, .x, sparse = FALSE)[, 1] | st_within(vect_id$geometry, .x, sparse = FALSE)[, 1] )
    })) %>%
    dplyr::select(-geometry)
  
  return(out)
  
}


################################### DATA #######################################
# Read and inspect data
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

vector_net <- read_sf('./spatial_data/VectornetMAPforMOODjan21.shp', 
                      crs = st_crs(nuts0)) %>%
  st_make_valid()

vector_net %<>% 
  st_transform(st_crs(vector_net))

vect_id <- as_tibble(vector_net[-1083,]) %>%
  rowid_to_column() %>%
  st_as_sf()

nuts_all <- bind_rows(nuts0,
                      nuts1,
                      nuts2)

metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')
clusters <- read_csv('./2025Jun24/europe_clusters/all_clusters.csv') %>%
  rename(tipnames = label)

seqid_cluster <- metadata_with_concat %>%
  left_join(clusters) %>%
  mutate(seq_id = coalesce(sequence_accession, sequence_isolate) %>%
           gsub('\\/', '_', .)) %>%
  select(seq_id, cluster)

################################### MAIN #######################################
# Main analysis or transformation steps

# For isolates with NUTS-0,1,2 spatial resolution, distribute proportionally 
# across all applicable polygons
list_uncertainty_all <- GetUncertaintyWeightings(metadata_with_concat, 
                                                 nuts_all) %>%
  left_join(seqid_cluster)

list_uncertainty_nflg <- GetUncertaintyWeightings(metadata_with_concat, 
                                                  nuts_all,
                                                  nflg_only = TRUE) %>%
  left_join(seqid_cluster)


# polygon number and how many assigned to each isolate (all)
uncertainty_weightings_all <- list_uncertainty_all %>%
  rowwise() %>%
  mutate(eurostat_polygon_n = length(eurostat_polygon)) %>%
  as_tibble() %>%
  mutate(weight_n = 1/eurostat_polygon_n) %>%
  unnest(eurostat_polygon) %>%
  select(rowid = eurostat_polygon,
         cluster,
         weight_n) %>%
  mutate(n = 1) %>%
  summarise(n = sum(n*weight_n), .by = c(rowid,cluster))


# polygon number and how many assigned to each isolate (nflg)
uncertainty_weightings_nflg <- list_uncertainty_nflg %>%
  rowwise() %>%
  mutate(eurostat_polygon_n = length(eurostat_polygon)) %>%
  as_tibble() %>%
  mutate(weight_n = 1/eurostat_polygon_n) %>%
  unnest(eurostat_polygon) %>%
  select(rowid = eurostat_polygon,
         cluster,
         weight_n) %>%
  mutate(n = 1) %>%
  summarise(n = sum(n*weight_n), .by = c(rowid,cluster))
            
# Summarise map metadata for nflg only
map_metadata_nflg <- metadata_with_concat %>%
  filter(generegion_nflg ==1) %>%
  left_join(clusters) %>%
  # format geo - extract coords from character string, assign separate cols and 
  # map to VectorNet polygons
  mutate(geocode_coords =  gsub('c\\(|\\)', '', geocode_coords)) %>%
  separate(geocode_coords, into = c("lat", "lon"), sep = ", ", convert = TRUE) %>%
  st_as_sf(coords = c('lat','lon'), crs = st_crs(vector_net), sf_column_name = 'coords') %>%
  st_join(.,
          vect_id, 
          join = st_within, 
          # left = FALSE,
          # largest = TRUE
  ) %>% 
  drop_na(nuts2_id) %>%
  
  # format
  count( rowid, cluster) %>%
  st_drop_geometry() 


# Summarise map metadata for all
map_metadata_all <- metadata_with_concat %>%
  left_join(clusters) %>%
  # format geo - extract coords from character string, assign separate cols and 
  # map to VectorNet polygons
  mutate(geocode_coords =  gsub('c\\(|\\)', '', geocode_coords)) %>%
  separate(geocode_coords, into = c("lat", "lon"), sep = ", ", convert = TRUE) %>%
  st_as_sf(coords = c('lat','lon'), crs = st_crs(vector_net), sf_column_name = 'coords') %>%
  st_join(.,
          vect_id, 
          join = st_within, 
          # left = FALSE,
          # largest = TRUE
  ) %>% 
  
  drop_na(nuts2_id) %>%
  
  # format
  count(rowid, cluster) %>%
  st_drop_geometry() 


# create map containing all polygons (nflg)
nflg_map_data <- expand_grid(#lineage = level_1_lineages,
  rowid = 1:1489,
  cluster = LETTERS[1:8],
  n = NA_integer_)  %>% 
  rows_update(map_metadata_nflg %>% drop_na(), by = c('rowid', 'cluster')) %>% 
  left_join(vect_id)  %>%
  left_join(uncertainty_weightings_nflg, by = c('rowid', 'cluster')) %>% 
  replace_na(list('n.x' = 0, 'n.y' = 0))  %>%
  mutate(n= n.x + n.y) %>%
  select(rowid, n, geometry, cluster) %>%
  mutate(n = if_else(n ==0, NA, n)) %>%
  st_as_sf()


# create map containing all polygons (all)
all_map_data <- expand_grid(#lineage = level_1_lineages,
  rowid = 1:1489,
  cluster = LETTERS[1:8],
  n = NA_integer_)  %>% 
  rows_update(map_metadata_all %>% drop_na(), by = c('rowid', 'cluster')) %>% 
  left_join(vect_id)  %>%
  left_join(uncertainty_weightings_nflg, by = c('rowid', 'cluster')) %>% 
  replace_na(list('n.x' = 0, 'n.y' = 0))  %>%
  mutate(n= n.x + n.y) %>%
  select(rowid, n, geometry, cluster) %>%
  mutate(n = if_else(n ==0, NA, n)) %>%
  st_as_sf()
# filter(CNTR_CODE != 'TR') 


plt_a <- nflg_map_data %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'grey80') + 
  #geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34.4,65), xlim = c(-11, 34.2), expand = FALSE) +
  #scale_fill_brewer(palette = 'RdPu')
  scale_fill_distiller(palette = 'RdPu', 
                       transform = 'log10', 
                       direction = 1, 
                       breaks = c(1*10**(-1), 1*10**0, 1*10**1, 1*10**2),
                       na.value="grey95") +
  scale_alpha_continuous(range = c(0.5, 1), na.value = 1) +
  facet_grid(cols = vars(cluster)) +
  theme_void() + 
  theme(legend.position = 'none')


plt_b <- all_map_data %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'grey80') + 
  #geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34.4,65), xlim = c(-11, 34.2), expand = FALSE) +
  #scale_fill_brewer(palette = 'RdPu')
  scale_fill_distiller(palette = 'RdPu', 
                       transform = 'log10', 
                       direction = 1, 
                       breaks = c(1*10**(-1), 1*10**0, 1*10**1, 1*10**2),
                       na.value="grey95") +
  scale_alpha_continuous(range = c(0.5, 1), na.value = 1) +
  facet_grid(cols = vars(cluster)) +
  theme_void() + 
  theme(legend.position = 'none')

################################### OUTPUT #####################################
# Save output files, plots, or results
cowplot::plot_grid(plt_a, plt_b, nrow = 2)


#################################### END #######################################
################################################################################