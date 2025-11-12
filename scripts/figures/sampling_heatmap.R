################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-10-21
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
library(lubridate)

################################### DATA #######################################
# Read and inspect data
metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv')

################################### MAIN #######################################
# Main analysis or transformation steps
my_theme <- theme_classic()+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


metadata %>%
  #mutate(is_europe = ifelse(is_europe == 1, 'Europe', 'Non-Europe')) %>%
  filter(is_europe == 1) %>%
  group_by(date_y, nuts0_id, is_europe) %>%
  summarise(n = n()) %>%
  ggplot() + 
  geom_tile(aes(x = date_y, y = nuts0_id, fill=n)) + 
  #scale_x_date('Collection Year') + 
  scale_y_discrete('Country') + 
  scale_fill_distiller(palette = 'OrRd', direction = 1) + 
  my_theme + 
  #facet_grid(rows = vars(is_europe), scales = 'free_y', space = 'free_y',switch="both") + 
  theme(legend.position = 'bottom', strip.placement = 'outside')
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################
