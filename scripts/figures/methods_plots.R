################################################################################
## Script Name:        Figure 3 Composite
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-12-01
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

nuts_id <- nuts0 %>%
  dplyr::select(geometry) %>%
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
  left_join(uncertainty_weightings_all, by = c('rowid', 'cluster')) %>% 
  replace_na(list('n.x' = 0, 'n.y' = 0))  %>%
  mutate(n= n.x + n.y) %>%
  select(rowid, n, geometry, cluster) %>%
  mutate(n = if_else(n ==0, NA, n)) %>%
  st_as_sf()
# filter(CNTR_CODE != 'TR') 

################################### OUTPUT #####################################
# Save output files, plots, or results
a <-  metadata_noconcat %>%
  filter(generegion_nflg ==1) %>%
  filter(is_europe == 1) %>%
  filter(drop_fli == 0) %>% #Problematic FLI sequences removed
  count(nuts0_id) %>%
  arrange(desc(n)) %>%
  mutate(nuts0_id = factor(nuts0_id, levels=unique(nuts0_id[order(n, decreasing = T)]), ordered = T)) %>%
  ggplot() +
  geom_bar(aes(x = nuts0_id, y = n), stat = 'identity') + 
  scale_x_discrete('Country') + 
  scale_y_continuous('N', expand = c(0,0), ) +
  theme_classic(base_size = 9)


# Plt - sequence length
b <- metadata_noconcat %>%
  ggplot()+
  #geom_rect(aes(xmin = 9000, xmax = Inf, ymin = 0, ymax = Inf), alpha = 0.007, fill = 'lightgrey')+
  geom_histogram(aes(x = sequence_length), binwidth = 150) +
  #geom_vline(xintercept  = 9000, linetype = 'dashed', linewidth = 1) + 
  geom_vline(xintercept  = 200, linetype = 'dashed', linewidth = 1) + 
  scale_x_continuous('Sequence Length', expand = c(0,0), breaks = seq(0,10000, by= 2000)) + 
  scale_y_continuous('Count', expand = c(0,0)) +
  scale_fill_brewer(palette = 'PuBu', na.value = 'black', 'Genomic Region') + 
  theme_classic(base_size = 9)


# Plt - alignment
c <- metadata_noconcat %>%
  filter(generegion_nflg == 0) %>%
  arrange(sequence_start) %>%
  dplyr::select(tipnames, sequence_start, sequence_end, nuts0_id) %>%
  pivot_longer(starts_with('sequence'), values_to = 'position', names_to = 'name') %>%
  dplyr::select(-name) %>%
  group_by(tipnames) %>%
  complete(position = seq(min(position), max(position))) %>%
  fill(nuts0_id, .direction = 'down') %>%
  fill(nuts0_id, .direction = 'down') %>%
  ungroup() %>%
  summarise(prop_coverage = n()/856, .by = position) %>%
  arrange(position) %>%
  ggplot() + 
  geom_ribbon(aes(x = position, ymax = prop_coverage, ymin = 0),
              colour = 'black',
              fill = 'blue',
              alpha = 0.5) +
  #geom_path(aes(x = position, y = prop_coverage), linewidth = 1) +
  scale_y_continuous('Proportion Coverage', expand = c(0,0)) + 
  scale_x_continuous('Alignment Position', expand= c(0.01, 0)) + 
  coord_cartesian(xlim = c(0,10000)) +
  theme_classic(base_size = 9)

top <- plot_grid(a, b, c, ncol = 3, align = 'v', axis = 'tb', labels = 'AUTO', label_size = 10, scale = 0.95)


#################################### END #######################################
################################################################################
#  clustersmethod, sampling heatmaps with/without partial sequences
p1 <- ggdraw() + draw_image("./2025Jun24/plots/usutu_methods.png", scale = 1)


plt_a <- nflg_map_data %>%
  filter(cluster != 'A') %>%
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

# facet these two
# decrease spatial scale
# Add heatmap and sequence length plot
plt_d <-  nflg_map_data %>%
  mutate(data = 'NFLG') %>%
  bind_rows(mutate(all_map_data, data = 'All'))  %>%
  mutate(data = factor(data, levels = c('NFLG', 'All'), ordered = T)) %>%
  filter(cluster %in% c('C', 'E', 'F')) %>%
  rowwise() %>%
  mutate(cluster = seq_along(LETTERS)[LETTERS == cluster] %>% as.roman()) %>%
  as_tibble() %>%
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
  facet_grid(cols = vars(as.character(cluster)),
             rows = vars(data),
             #labeller = as_labeller(str_to_upper),
             switch = 'y') +
  theme_void() + 
  theme(legend.position = 'none')

lhs <- align_plots(a, p1, align = 'h', axis = 'l')
top <- plot_grid(lhs[[1]], b, c, ncol = 3, align = 'v', axis = 'tb', labels = 'AUTO', label_size = 10, scale = 0.95)

plot_grid(top, lhs[[2]], nrow = 2, labels = c('','D'), label_size = 10)

ggsave('./2025Jun24/plots/figure3.jpeg',
       dpi = 360,
       height = 15,
       width = 20,
       units = 'cm')