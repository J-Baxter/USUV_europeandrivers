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
temp <- read_sf('~/Downloads/NUTS_RG_10M_2021_4326.shp/NUTS_RG_10M_2021_4326.shp')

nongisco_shapefile <- ne_countries(scale = 10, country = c('bosnia and herzegovina', 'kosovo', 'andorra'), returnclass = "sf") %>%
  dplyr::select(name_en,
         iso_a2_eh,
         geometry) %>%
  rename(CNTR_CODE = iso_a2_eh,
         NAME_LATN = name_en) %>%
  mutate(NUTS_ID = CNTR_CODE)

ggplot(balkan_noneu) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE) + 
  theme_void()

regions %>%
  select(iso_a2, geometry)


italy_3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "10", #1:10million
  nuts_level = "3",
  country = 'italy') 

italy_3 %>%
  mutate(test = case_when(grepl('Bologna', NUTS_NAME) ~ '1', .default = 0))

ggplot(italy_3) + 
  geom_sf(aes(fill = grepl('Bologna', NUTS_NAME))) + 
  #coord_sf(ylim = c(34,60), xlim = c(-12, 45), expand = FALSE) + 
  theme_void()

nuts2 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "10", #1:10million
  nuts_level = "2")

nuts0 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "10", #1:10million
  nuts_level = "0")

all_europe_sf <- bind_rows(nuts2, nongisco_shapefile)

ggplot(all_europe_sf) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE) + 
  theme_void()



# Import data
metadata <- read_csv('./data/USUV_metadata_noFLI_2024Aug28.csv')
lu_spatialdata <- readRDS('./spatial_data/Predictors_tep_20231014.Rdata')

#test raster
test_df <- lu_spatialdata %>% 
  dplyr::select(c(lon, lat, Deciduous_broad))

# metadat
test_metadata <- metadata %>%
  drop_na(c('collection_subdiv1lat', 'collection_subdiv1long')) %>%
  st_as_sf(coords = c('collection_subdiv1long', 'collection_subdiv1lat'), crs = 4326) %>%
  st_join(nuts2, join = st_within, largest = TRUE,  left = FALSE) %>%
  dplyr::select(-c(
    starts_with('collection_subdiv1'),
    starts_with('collection_subdiv2')
  )) %>%
  dplyr::select(-geometry) %>%
  as_tibble()


plt1 <- expand_grid(sequence_generegion = unique(test_metadata$sequence_generegion),
            NUTS_ID = all_europe_sf$NUTS_ID,
            n = NA_integer_) %>%
  rows_patch(test_metadata %>%
                summarise(n = n(), .by = c(NUTS_ID, sequence_generegion)),
              by = c('NUTS_ID', 'sequence_generegion')) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', transform = 'log10', direction = -1, , na.value="lightgrey", breaks = c(0, 1, 5, 10, 20, 50, 100)) +
  facet_grid(cols = vars(sequence_generegion)) +
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
expand_grid(sequence_generegion = unique(test_metadata$sequence_generegion),
            NUTS_ID = all_europe_sf$NUTS_ID,
            n = NA_integer_,
            date_y = ) %>%
  rows_patch(test_metadata %>%
               summarise(n = n(), .by = c(NUTS_ID, sequence_generegion, date_y)),
             by = c('NUTS_ID', 'sequence_generegion')) %>%
  left_join(all_europe_sf, .) %>%
  filter(CNTR_CODE != 'TR') %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', transform = 'log10', direction = -1, , na.value="lightgrey", breaks = c(0, 1, 5, 10, 20, 50, 100)) +
  facet_grid(cols = vars(sequence_generegion)) +
  theme_void() + 
  theme(legend.position = 'none')
ggplot(nuts2) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE) + 
  geom_raster(aes(y = lat, x = lon, fill = Deciduous_broad), data = test_df)


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

