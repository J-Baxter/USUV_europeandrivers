################################################################################
## Script Name:        Plot Skygrids
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-08-21
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
library(treeio)
library(beastio)

################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]

skygridpaths <- sapply(dirs, 
                       list.files,
                       pattern = "\\.skygrid$",
                       full.names = TRUE, 
                       simplify = F) %>%
  Filter(length,.)

skygrid_logs <- lapply(skygridpaths,
                       read_delim, 
                       delim ='\t', 
                       skip = 1) %>%
  bind_rows(., .id = 'clade') 

test_grid <- read_delim('./2025Jun24/europe_clusters/NFLG_V/Untitled', delim ='\t', skip = 1) 
################################### MAIN #######################################
# Main analysis or transformation steps

skygrid_logs %>%
  mutate(clade = gsub('\\.\\/2025Jun24\\/europe_clusters\\/', '', clade)) %>%
  separate_wider_delim(clade, delim  = '_', names = c('sequence', 'clade')) %>%
  ggplot() + 
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper), alpha = 0.5) + 
  geom_line(aes(x = date, y=median)) + 
  scale_x_date(expand = c(0,0)) + 
  scale_y_log10(expand = c(0,0)) + 
  theme_classic() + 
  facet_grid(cols = vars(sequence), rows = vars(clade), scales = 'free')


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################
# Deprecated
#skygrid_b <- read_delim('./2024Oct20/alignments/subset_alignments/run/B_skygrid.txt') %>%
#  separate_wider_delim(5, '\t', names = c('mean', 'median', 'upper', 'lower')) %>% 
#  .[-1,] %>%
#  rename(time = ...1) %>%
#  mutate(across(c(1,5,6,7,8), .fns = ~as.numeric(.x)))

#ggplot(skygrid_b) +
#  geom_line(aes(x = time, y = median), colour = 'brown') + 
#  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = 'brown', alpha = 0.1)+
#  scale_y_log10()+
#  theme_minimal()