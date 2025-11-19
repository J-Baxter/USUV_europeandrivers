pasture_1kmsq <- rast('./2025Jun24/raster_data/pastures.tiff') %>%
  as.int() %>%
  aggregate(., fact=10, fun = 'mean')

moors_and_heathland_1kmsq <- rast('./2025Jun24/raster_data/moors_and_heathland.tiff') %>%
  as.int() %>%
  aggregate(., fact=10, fun = 'mean')

discontinuous_ubran_1kmsq <- rast('./2025Jun24/raster_data/discontinuous_urban_fabric.tiff') %>%
  as.int() %>%
  aggregate(., fact=10, fun = 'mean')

coniferous_1kmsq <- rast('./2025Jun24/raster_data/coniferous_forest.tiff') %>%
  as.int() %>%
  aggregate(., fact=10, fun = 'mean')



plt_pasture <- ggplot() +
  geom_spatraster(data = pasture_1kmsq,  maxcell = 1e+06) +
  scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'grey99')+
  theme_void()+
  theme(legend.position = 'None',
        legend.title = element_blank())

plt_moors <- ggplot() +
  geom_spatraster(data = moors_and_heathland_1kmsq,  maxcell = 1e+06) +
  scale_fill_distiller(palette = 'Reds', direction = 1, na.value = 'grey99')+
  theme_void()+
  theme(legend.position = 'None',
        legend.title = element_blank())




mdt <- get_elev_raster(discontinuous_ubran_1kmsq, z = 5)
mdt <- project(rast(mdt), crs(vector_net))


# convert to terra and mask area of interest
ext(mdt) <- ext(discontinuous_ubran_1kmsq)
mdt <- resample(mdt, discontinuous_ubran_1kmsq)
mdt_masked <- mask(mdt, discontinuous_ubran_1kmsq)

# reproject

sl <- terrain(mdt_masked, "slope", unit = "radians")
asp <- terrain(mdt_masked, "aspect", unit = "radians")

hill_single <- shade(sl, asp,
                     angle = 45,
                     direction = 300
)

plot(hill_single)

test <- c(hill_single, discontinuous_ubran_1kmsq)

plt_discontinuous_urban <- ggplot() +
  geom_spatraster(data = hill_single,  maxcell = 1e+06) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_spatraster(data = discontinuous_ubran_1kmsq,  maxcell = 1e+0) +
  scale_fill_distiller(palette = 'Blues', direction = 1, na.value = '#FF000000')+
  theme_void()+
  theme(legend.position = 'None',
        legend.title = element_blank())

hill_plot <- ggplot() +
  geom_spatraster(data = hill_single,  maxcell = 1e+07, aes(fill = hillshade)) +
  scale_fill_distiller(palette = 'Greys', direction = 1, na.value = 'grey97') +
  theme_void()+
  theme(legend.position = 'None',
        legend.title = element_blank())


hill_plot + 
  new_scale_fill() +
  geom_spatraster(data = pasture_1kmsq, maxcell = 1e+07, alpha = 0.75) +
  scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'grey97') +
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) 

hill_plot + 
  geom_spatraster(data = coniferous_1kmsq, maxcell = Inf, alpha = 0.4) +
  scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'grey97') 




plt_coniferous <- ggplot() +
  geom_spatraster(data = coniferous_1kmsq,  maxcell = 1e+06) +
  scale_fill_distiller(palette = 'Greys', direction = 1, na.value = 'grey99')+
  theme_void()+
  theme(legend.position = 'None',
        legend.title = element_blank())

cowplot::plot_grid(plt_pasture, plt_moors, plt_discontinuous_urban, plt_coniferous,
                   labels = c('Pasture', 'Moorland ', 'Urban', 'Coniferous Forest'),
                   nrow = 1)
