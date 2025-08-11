

culex_abundance_files <- list.files('./culex_models/statistical_models/Abundace_model_predictions',
                                    full.names = T)

culex_cs_files <- list.files('./culex_models/statistical_models/Abundance_plus_CitizenScience_model_predictions',
                                    full.names = T)


culex_abundance <- lapply(culex_abundance_files, read_csv) %>%
  setNames(gsub('\\.\\/culex_models\\/statistical_models\\/Abundace_model_predictions\\/culex_|\\.csv', '', culex_abundance_files)) %>%
  bind_rows(., .id = 'mm_yyyy') %>%
  separate_wider_delim(., mm_yyyy, '_', names = c('month', 'year')) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) 

culex_cs <- lapply(culex_cs_files, read_csv) %>%
  setNames(gsub('\\.\\/culex_models\\/statistical_models\\/Abundance_plus_CitizenScience_model_predictions\\/culex_|\\.csv|_ma', '', culex_cs_files)) %>%
  bind_rows(., .id = 'mm_yyyy') %>%
  separate_wider_delim(., mm_yyyy, '_', names = c('month', 'year')) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) 

culex_abundance_aggregate <- aggregate(culex_abundance %>% 
                                         dplyr::select(prediction,geometry),
                                       vect_id, median)

culex_cs_aggregate <- aggregate(culex_cs %>% 
                                  dplyr::select(prediction,geometry),
                                vect_id, mean)



ggplot(culex_abundance_by_vectornet[[1]] , aes(colour = prediction)) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_colour_distiller(palette = 'RdYlGn', 
                       #transform = 'log10', 
                       direction = -1, , 
                       na.value="lightgrey") +
 # facet_grid(rows = vars(year), cols = vars(month)) + 
  theme_void() 



ggplot(culex_abundance_aggregate,
       aes(fill = prediction)) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdYlGn', 
                         #transform = 'log10', 
                         direction = -1, , 
                         na.value="lightgrey") +
  #facet_grid(rows = vars(year), cols = vars(month)) + 
  theme_void() 
