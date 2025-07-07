################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-06-04
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
#library(eurostat)


################################### DATA #######################################
# Read and inspect data

nongisco_shapefile <- ne_countries(scale = 10, country = c('bosnia and herzegovina', 'kosovo', 'andorra'), returnclass = "sf") %>%
  dplyr::select(name_en,
                iso_a2_eh,
                geometry) %>%
  rename(CNTR_CODE = iso_a2_eh,
         NAME_LATN = name_en) %>%
  mutate(NUTS_ID = CNTR_CODE)

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
  nuts_level = "3") %>%
  ggplot() + geom_sf() +  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE)


nuts3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "3")

all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)

vector_net <- read_sf('./spatial_data/VectornetMAPforMOODjan21.shp', crs = st_crs(nuts0)) %>%
  st_make_valid() %>%
  st_transform(st_crs(vector_net))

metadata <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')

################################### MAIN #######################################
# Plot Maps of Europe stratified by selected lineage
level_1_lineages <- lineage_tbl$`GRI Lineage Level 1` %>% unique() %>% .[!is.na(.)]


vect_id <- as_tibble(vector_net[-1083,]) %>%
  rowid_to_column() %>%
  st_as_sf()

map_metadata <- metadata %>%
  
  # format geo
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
  left_join(lineage_tbl, by = join_by('tipnames' == 'label')) %>%
  drop_na(`GRI Lineage Level 1`) %>%
  count(`GRI Lineage Level 1`, rowid) %>%
  rename(lineage = `GRI Lineage Level 1`) %>%
  st_drop_geometry()
  
  
map_data <- expand_grid(lineage = level_1_lineages,
            rowid = vect_id$rowid,
            n = NA_integer_) %>%
  rows_patch(map_metadata,
             unmatched = 'ignore',
             by = c('rowid', 'lineage')) %>%
  left_join(vect_id) %>%
  filter(grepl('^A', lineage)) 
 # filter(CNTR_CODE != 'TR') 


map_data %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  #geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', 
                       #transform = 'log10', 
                       direction = -1, , 
                       na.value="lightgrey", 
                       breaks = c(0, 1, 5, 10, 20, 50, 100)) +
  facet_wrap(~lineage) +
  theme_void() + 
  theme(legend.position = 'none')




mdt <- get_elev_raster(vector_net, z = 5)

# convert to terra and mask area of interest
mdt <- rast(mdt) |>
  mask(vect(vector_net))

# reproject
mdt <- project(mdt, crs(vector_net))

sl <- terrain(mdt, "slope", unit = "radians")
asp <- terrain(mdt, "aspect", unit = "radians")

hill_single <- shade(sl, asp,
                     angle = 45,
                     direction = 300,
                     normalize = TRUE
)

plot(hill_single)

hilldf_single <- as.data.frame(hill_single, xy = TRUE)
mdtdf <- as.data.frame(mdt, xy = TRUE)
names(mdtdf)[3] <- "alt"

ggplot() +
  geom_raster(
    data = hilldf_single,
    aes(x, y, fill = hillshade),
    show.legend = FALSE
  ) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(
    data = mdtdf,
    aes(x, y, fill = alt),
    alpha = .7
  ) +
  scale_fill_hypso_tint_c(breaks = c(
    180, 250, 500, 1000,
    1500, 2000, 2500,
    3000, 3500, 4000
  )) +
 
  guides(fill = guide_colorsteps(
    barwidth = 20,
    barheight = .5,
    title.position = "right"
  )) +
  labs(fill = "m") +
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  theme_void() + 
  theme(legend.position = 'none')




# reproject vect
suiz <- st_transform(suiz, st_crs(suiz_lakes))
plot(mdt)

vector_net <- spdep::poly2nb(st_make_valid(vector_net)) %>%
  st_transform(st_crs(vector_net))



  map %>%
  ggplot() + 
  geom_sf() + 
  geom_sf(data = nuts0, colour = 'black', linewidth = 0.2, fill = NA) +
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  theme_void()

metadata %>%
  mutate(geocode_coords =  gsub('c\\(|\\)', '', geocode_coords)) %>%
  separate(geocode_coords, into = c("lat", "lon"), sep = ", ", convert = TRUE) %>%
  st_as_sf(coords = c('lat','lon'), crs = st_crs(vector_net), sf_column_name = 'geometry') %>% 
  st_join(.,
          vector_net[-1083,], 
          join = st_within, 
          # left = FALSE,
         # largest = TRUE
          ) %>% 
  mutate(geometry = if_else(is.na(nuts3_id), NA, geometry)) %>%
  count()



################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################
# paper figures



# Initialise map of europe subdivisions (NUTS level 2)
# Data obtained from Eurostat GISCO (the Geographic Information System of the COmmission) via API 
# package giscoR
#temp <- read_sf('~/Downloads/NUTS_RG_10M_2021_4326.shp/NUTS_RG_10M_2021_4326.shp')



#ggplot(balkan_noneu) + 
  #geom_sf() + 
 # coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE) + 
 # theme_void()

#regions %>%
  #select(iso_a2, geometry)


#italy_3 <- gisco_get_nuts(
 # year = "2021",
 # epsg = "4326", #WGS84 projection
 # resolution = "10", #1:10million
 # nuts_level = "3",
 # country = 'italy') 

#italy_3 %>%
 # mutate(test = case_when(grepl('Bologna', NUTS_NAME) ~ '1', .default = 0))

#ggplot(italy_3) + 
 # geom_sf(aes(fill = grepl('Bologna', NUTS_NAME))) + 
  #coord_sf(ylim = c(34,60), xlim = c(-12, 45), expand = FALSE) + 
  #theme_void()




# Import data
lu_spatialdata <- readRDS('./spatial_data/Predictors_tep_20231014.Rdata')

#test raster
test_df <- lu_spatialdata %>% 
  dplyr::select(c(lon, lat, Deciduous_broad))

# metadat
#test_metadata <- metadata %>%
  #drop_na(c('collection_subdiv1lat', 'collection_subdiv1long')) %>%
  #st_as_sf(coords = 'geocode_coords', crs = 4326) %>%
  #st_join(nuts2, join = st_within, largest = TRUE,  left = FALSE) %>%
  #dplyr::select(-c(
    #starts_with('collection_subdiv1'),
    #starts_with('collection_subdiv2')
  #)) %>%
  #dplyr::select(-geometry) %>%
  #as_tibble()

all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)

plt1 <- expand_grid(gene_region = c('nflg', 
                                    'env',
                                    'NS5'),
                    NUTS_ID = all_europe_sf$NUTS_ID,
                    n = NA_integer_) %>%
  rows_patch(metadata %>%
               rename(NUTS_ID = nuts1_id) %>%
               pivot_longer(starts_with('generegion'), names_to = 'gene_region') %>%
               mutate(gene_region = gsub('^generegion_', '', gene_region)) %>%
               summarise(n = sum(value), .by = c(NUTS_ID, gene_region)),
             unmatched = 'ignore',
             by = c('NUTS_ID', 'gene_region')) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', 
                       #transform = 'log10', 
                       direction = -1, , 
                       na.value="lightgrey", 
                       breaks = c(0, 1, 5, 10, 20, 50, 100)) +
  facet_grid(cols = vars(gene_region)) +
  theme_void() + 
  theme(legend.position = 'none')
  
plt2 <- test_metadata %>%
  ggplot() + 
  geom_histogram(aes(x = date_ymd), binwidth = 365, fill = '#7a0177') + 
  scale_x_date(expand = c(0.02,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + 
  facet_grid(cols = vars(sequence_generegion)) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside")


cowplot::plot_grid(plt1, plt2, nrow = 2, align = 'hv')


# alternative plot (facet by genomic region and year)
all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)

expand_grid(gene_region = c('nflg', 
                            'env_1003_1491',
                            #'NS5_9100_9600', 
                            #'NS5_8968_9264',
                            'NS5_9000_9600',
                            'NS5_10042_10312'),
            NUTS_ID = all_europe_sf$NUTS_ID,
            year_group = c("<2004-2010", "2011-2015", "2016-2020", '2021-2024'),
            n = NA_integer_) %>%
  rows_patch(metadata %>%
               rename(NUTS_ID = nuts1_id) %>%
               pivot_longer(starts_with('generegion'), names_to = 'gene_region') %>%
               mutate(gene_region = gsub('^generegion_', '', gene_region)) %>%
               mutate(year_group = cut(as.numeric(date_y), 
                                       breaks = c(0, 2010, 2015, 2020, 2024), 
                                       labels = c("<2004-2010", "2011-2015", "2016-2020", '2021-2024'))) %>%
               summarise(n = sum(value), .by = c(NUTS_ID, gene_region, year_group)),
             unmatched = 'ignore',
             by = c('NUTS_ID', 'gene_region', 'year_group')) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white', linewidth = 0.1) + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.2, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', transform = 'log10', direction = -1, , na.value="lightgrey", breaks = c(0, 1, 5, 10, 20, 50, 100)) +
  facet_grid(rows = vars(gene_region),
             cols = vars(year_group)) +
  theme_void() + 
  theme(legend.position = 'none')


#year only
all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)
expand_grid( NUTS_ID = all_europe_sf$NUTS_ID,
            year_group = c("<2004-2010", "2011-2015", "2016-2020", '2021-2024'),
            n = NA_integer_) %>%
  rows_patch(metadata %>%
               rename(NUTS_ID = nuts1_id) %>%
               mutate(year_group = cut(as.numeric(date_y), 
                                       breaks = c(0, 2010, 2015, 2020, 2024), 
                                       labels = c("<2004-2010", "2011-2015", "2016-2020", '2021-2024'))) %>%
               summarise(n = n(), .by = c(NUTS_ID, year_group)),
             unmatched = 'ignore',
             by = c('NUTS_ID', 'year_group')) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white', linewidth = 0.1) + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.2, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_fermenter(palette = 'RdPu',
                       direction = -1, 
                       transform = 'log10',
                       na.value="lightgrey", 
                       breaks = c(1, 5, 10, 20, 50, 100)) +
  facet_grid(cols = vars(year_group)) +
  theme_void(base_size = 20) + 
  theme(legend.position = 'bottom') 


coordinates(test_df) <- ~ lon + lat
r <- raster(extent(test_df), resolution = 0.008333333)
crs(r) <- CRS("+proj=longlat +datum=WGS84")

raster_test <- rasterize(test_df, r, field = 'Deciduous_broad', fun = mean)

chicken <- raster('~/Downloads/Chickens_Mekong.asc/Chickens_Mekong.asc')
plot(chicken)

chicken <- raster('~/Downloads/GLW4-2020.D-DA.CHK.tif')
plot(chicken)

cattle <- raster('~/Downloads/GLW4-2020.D-DA.CTL.tif')
plot(log10(cattle))

######################## 30th October
clusters <- read_csv('./data/clusters_2024Oct29.csv')
concat_metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv') %>%
  left_join(clusters)
ml_concat <- read.tree('./2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_withconcatenated.fasta.contree')

# first root ml tree by south african sequences
rooted_concat <- root(ml_concat, outgroup =  "EU074021|NA|NA|ZA|1958" , resolve.root = TRUE)


# Determine the MRCA of clusters, as determined by NFLG sequences
cluster_mrcas <- rooted_concat %>%
  as.treedata() %>%
  left_join(concat_metadata,
            by = join_by(label == tipnames)) %>%
  as_tibble() %>%
  drop_na(cluster) %>%
  group_split(cluster) %>%
  lapply(., pull, 'node') %>%
  lapply(., getMRCA, phy = rooted_concat)


cluster_offspring <- lapply(cluster_mrcas, offspring, .data =  rooted_concat, tiponly = TRUE) %>%
  lapply(., function(x) rooted_concat$tip.label[x]) %>% 
  lapply(., as_tibble) %>%
  setNames(LETTERS[1:5]) %>%
  bind_rows(., .id = 'cluster') %>%
  rename(tipnames = value)

cluster_metadata <- concat_metadata %>%
  rows_patch(., cluster_offspring, by='tipnames') 

cluster_metadata_split <- cluster_metadata %>%
  group_split(cluster)


# Requires cluster designation for partial genomes
nlfg_mapdata <- expand_grid(gene_region = 'nflg',
                            cluster = LETTERS[1:5],
                            NUTS_ID = all_europe_sf$NUTS_ID,
                            n = NA_integer_) %>%
  rows_patch(metadata %>%
               left_join(clusters) %>%
               rename(NUTS_ID = nuts1_id) %>%
               filter(generegion_nflg == 1) %>%
               pivot_longer(starts_with('generegion'), names_to = 'gene_region') %>%
               mutate(gene_region = gsub('^generegion_', '', gene_region)) %>%
               summarise(n = sum(value), .by = c(NUTS_ID, gene_region, cluster)),
             unmatched = 'ignore',
             by = c('NUTS_ID', 'gene_region', 'cluster')) %>%
  mutate(row = 'nflg')

all_sequences_mapdata <- expand_grid(cluster = LETTERS[1:5],
                                     NUTS_ID = all_europe_sf$NUTS_ID,
                                     n = NA_integer_) %>%
  rows_patch(cluster_metadata %>%
               rename(NUTS_ID = nuts1_id) %>%
               summarise(n = n(), .by = c(NUTS_ID, cluster)),
             unmatched = 'ignore',
             by = c('NUTS_ID', 'cluster')) %>%
  mutate(row = 'all_sequences')

bind_rows(nlfg_mapdata,
          all_sequences_mapdata) %>%
  mutate(row = factor(row, levels = c('nflg', 'all_sequences'))) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_fermenter(palette = 'RdPu',
                       direction = -1, 
                       #transform = 'log10',
                       na.value="lightgrey", 
                       breaks = c(1, 5, 10, 20, 50, 100)) +
  facet_grid(cols = vars(cluster), rows= vars(row),
             labeller = labeller(cluster = c('A' = 'A (EU1+2)',
                                             'B' = 'B (EU3)',
                                             'C' = 'C (AF3)',
                                             'D' = 'D (AF2)',
                                             'E' = 'E (AF1)'),
                                 row = c('all_sequences' = 'All sequences',
                                         'nflg' = 'NFLG')), switch = 'y') +
  theme_void(base_size = 18) 
