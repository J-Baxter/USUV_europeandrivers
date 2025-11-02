################################################################################
## Script Name:        Global Map of Usutu Cases & Sequences
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-10-31
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
library(countrycode)


################################### DATA #######################################
# Read and inspect data
map <- ne_countries(scale = "medium", returnclass = "sf")

metadata <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')

usuv_pos_countries <- read_csv('./data/usuv_pos_countries.csv')

################################### MAIN #######################################
# Main analysis ors transformation steps
seqs_per_country <- metadata %>%
  count(nuts0_id) %>% 
  left_join(codelist %>% select(iso3c, nuts0_id = eurostat)) 

map_data <- map %>%
  left_join(usuv_pos_countries, 
            by = join_by(iso_a3_eh == iso3c)) %>%
  left_join(seqs_per_country, 
            by = join_by(iso_a3_eh == iso3c))
  

usuv_global_map <-  ggplot() +
  geom_sf_pattern(data = map_data %>%
                    mutate(usuv_positive = as.character(usuv_positive)) %>%
                    drop_na(usuv_positive),
    aes(pattern = usuv_positive),
                  na.rm = TRUE,
                  pattern_density = 0.01,
                  pattern_spacing = 0.01,
                  pattern_colour = 'red') +
  scale_pattern_manual(values = "stripe", labels = 'USUV Detected') +
  geom_sf(data = map_data,
          aes(fill =n),
          inherit.aes = F) +
  scale_fill_distiller(palette = 'Reds', direction = 1, na.value = '#FF000000') + 
  theme_void() +
  coord_sf(ylim = c(-38, 70),
           xlim = c(-20, 50),
           expand = TRUE) + 
  guides(pattern = guide_legend(NULL, position = 'inside', order = 1),
         fill = guide_colourbar('Sequences (n)', position = 'inside', order = 2)) + 
  theme(panel.background = element_rect(fill = "grey99",
                                        colour = '#FF000000'),
        legend.position.inside = c(0.15,0.15))


ggsave('./2025Jun24/plots/usuv_global_map.jpeg', height = 30, width = 20, units = 'cm', dpi =360)

  
################################### OUTPUT #####################################
# Save output files, plots, or results
#################################### END #######################################
################################################################################