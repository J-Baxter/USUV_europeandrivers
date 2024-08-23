# paper figures
library(sf)
library(rnaturalearth)
library(tidyverse)


col <- base::unique(data$colour[!is.na(data$colour)])

# Load base map
map <- ne_countries(returnclass = "sf", scale = 'medium') %>%
  left_join(metadata, by = join_by(adm0_iso), multiple = "all")

# Load subdisivion data
regions <- ne_states(country = c('United Kingdom', 'France', 'Portugal', 
                                 'Spain', 'Italy', 'Netherlands', 'Germany', 
                                 'Belgium','Sweden', 'Norway', 'Ireland', 
                                 'Poland', 'Czech Republic', 'Austria',
                                 'Slovenia', 'Switzerland', 'Denmark',
                                 'Hungary','Romania', 'Greece'), returnclass = "sf") %>%
  
# Plot maplibr
# Coordinates for Europe (Excl Iceland) ylim = c(35,72), xlim = c(-12, 45)

plot <- ggplot() +
  geom_sf(data = map,linewidth = 0.1) + 
  geom_sf(data = regions,linewidth = 0.1) + 
  coord_sf(ylim = c(35,72), xlim = c(-12, 45), expand = FALSE) +
  geom_point(data = metadata, aes(x = collection_subdiv1long, collection_subdiv1lat))
  my_theme + 
  scale_fill_gradientn(colors = create_sequential_palette(col, 11)[3:10] ,
                       limits = c(0,300),
                       breaks = seq(0,300,by=75),
                       na.value=create_sequential_palette(col, 10)[2]) +
  facet_grid(cols = vars(flu_season),
             drop = F)+
  my_theme +
  theme(legend.position = 'right',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  guides(fill = guide_coloursteps(barwidth = 0.76, barheight = 4))