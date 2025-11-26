vector_net <- read_sf('./spatial_data/VectornetMAPforMOODjan21.shp', crs = 4326) %>%
  st_make_valid() %>%
  st_transform(3035)

mdt <- get_elev_raster(pastures, src = 'gl3')

ebba_shapefile <- read_sf('./2025Jun24/raster_data/ebba/ebba2_etrs10x10_v1/ebba2_etrs10x10_v1.shp') %>%
  st_transform(4326)

nuts0 <- gisco_get_nuts(
  year = "2021",
  epsg = "3035", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "0") 

pastures <- rast('./2025Jun24/raster_data/pastures.tiff') %>%
  as.int()

pastures_agg <- aggregate(pastures, fact=2, fun = 'mean')

vals <- terra::extract(pastures_agg, 
                       ebba_shapefile,
                       fun = mean, 
                       na.rm = TRUE)


ebba_shapefile$pastures <- vals[,2]  
ebba_shapefile %>%
  ggplot(., aes(fill = pastures)) +
  geom_sf(lwd = 0, colour = '#ff000000') +
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  theme_void()

ebba_csv <- read_delim('./2025Jun24/raster_data/ebba/ebba2_data_model_10km.csv', delim = ';')

t_merula <- ebba_shapefile %>%
  left_join(ebba_csv %>% filter(birdlife_scientific_name == birds[[6]])) 

for (i in 1:length(birds)){
  tmp <- ebba_shapefile %>%
    left_join(ebba_csv %>% filter(birdlife_scientific_name == birds[[i]])) 
  
  r_temp <- st_rasterize(tmp %>% dplyr::select(`Probability of occurrence`, geometry))
  
  filename <- paste0('./2025Jun24/raster_data/ebba/', str_to_lower(birds[[i]]) %>% str_replace(' ', '_'), '.tiff')

  write_stars(r_temp, filename)
}

# rasterize based on geometry and a column named "value". Change the name of this column if necessary

# export as tiff

birds <- unique(ebba_csv$birdlife_scientific_name)
ebba_shapefile %>%
  
  #  %>%
  ggplot( vals) +
  geom_sf(colour = '#ff000000') + 
  #scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'grey99') + 
  theme_void()

ebba_csv %>% distinct(birdlife_scientific_name)


blackbird_raster <-  ebba_shapefile %>%
  left_join(ebba_csv %>% filter(birdlife_scientific_name == 'Turdus merula')) 
  




# Raster (SpatRaster)

# Polygons (SpatVector)
v <- vect('./2025Jun24/raster_data/ebba/ebba2_etrs10x10_v1/ebba2_etrs10x10_v1.shp')





head(extract(r, v, na.rm = TRUE))
v$avg <- extract(r, v, mean, na.rm = TRUE)$elevation


pastures_test <- pastures %>% 
  st_as_sf()


ggplot() +
  geom_spatraster(data = blackbird_raster,  maxcell = 1e+06) +
  scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'grey99')+
  theme_void()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())


ebba_shapefile_vect <- vect(ebba_shapefile)

ex <- extract(pastures, ebba_shapefile_vect, fun=mean, na.rm=TRUE)
ebba_shapefile$mean_val <- ex$lyr.1




ebba_shapefile %>%
  left_join(ebba_csv %>% filter(birdlife_scientific_name == birds[[6]])) %>%
  vect()