# paper figures
library(terra)
library(sf)
library(rnaturalearth)
library(giscoR)
#library(eurostat)
library(tidyverse)


# Initialise map of europe subdivisions (NUTS level 2)
# Data obtained from Eurostat GISCO (the Geographic Information System of the COmmission) via API 
# package giscoR
#temp <- read_sf('~/Downloads/NUTS_RG_10M_2021_4326.shp/NUTS_RG_10M_2021_4326.shp')

nongisco_shapefile <- ne_countries(scale = 10, country = c('bosnia and herzegovina', 'kosovo', 'andorra'), returnclass = "sf") %>%
  dplyr::select(name_en,
         iso_a2_eh,
         geometry) %>%
  rename(CNTR_CODE = iso_a2_eh,
         NAME_LATN = name_en) %>%
  mutate(NUTS_ID = CNTR_CODE)

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

nuts3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "3")

all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)

# Import data
metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')
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
                                    'env_1003_1491',
                                    'NS5_9100_9600', 
                                    'NS5_8968_9264',
                                    'NS5_10042_10312'),
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

# NFLG only
expand_grid(gene_region = 'nflg',
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
  facet_grid(cols = vars(cluster)) +
  theme_void() 

# Requires cluster designation for partial genomes
expand_grid(gene_region = c('nflg', 
                            'env_1003_1491',
                            'NS5_9100_9600', 
                            'NS5_8968_9264',
                            'NS5_10042_10312'),
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
  facet_grid(cols = vars(cluster)) +
  theme_void() 
