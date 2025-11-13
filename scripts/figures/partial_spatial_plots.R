################################################################################
## Script Name:       Plot spatial statistics and phylogeography 
## Purpose:           Plot phylogeography on a map, alongside spatial statistics
##                    including wavefront distances, diffusion coefficients and 
##.                   isolation-by-distance (IBD) signal
## Author:             James Baxter
## Date Created:       2025-11-12
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
# Load required libraries
library(tidyverse)
library(magrittr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(diagram)
library(lubridate)
library(treeio)
library(scales)


ReadSpatialStatistics <- function(filepath){
  out <-  read_delim(filepath,
                     delim = '\t', 
                     name_repair = function(nm) {
                       nm[2:length(nm)] <- paste0(seq_len(length(nm) - 1))
                       nm}) %>%
    pivot_longer(-1, names_to = 'iter', values_to = 'draw') 
  
  return(out)
  
}

################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]

wavefront_files <- sapply(dirs, list.files, pattern = 'wavefront_distances.txt',
                          full.names = TRUE, 
                          simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()

diffusioncoef_files <- sapply(dirs, list.files, pattern = 'tree_weighted_diffusion_coefficient.txt',
                          full.names = TRUE, 
                          simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()

wavefront_tbl <- lapply(wavefront_files, ReadSpatialStatistics) %>%
  setNames(wavefront_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(metric = gsub('USUV_|s.txt', '', metric))

diffusioncoef_tbl <- lapply(diffusioncoef_files, ReadSpatialStatistics) %>%
  setNames(diffusioncoef_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(metric = gsub('USUV_tree_|\\.txt', '', metric))
  
wavefront_tbl %>%
  ggplot() + 
  geom_path(aes(x = time, y = draw, colour = clade, group = iter)) + 
  facet_grid(rows = vars(metric), cols = vars(clade))


diffusioncoef_tbl %>%
  ggplot() + 
  geom_path(aes(x = time, y = draw, colour = clade, group = iter)) + 
  facet_grid(rows = vars(clade), scales = 'free_y')


################################### MAIN #######################################
# Main analysis or transformation steps
# Plots
# 1. Spatial Wavefront Distance
# 2. Patristic Wavefront Distance
# 3. Weighted Diffusion Coefficient
# 4. Diffusion Coefficient Variation 
# 5. Isolate by distance
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################