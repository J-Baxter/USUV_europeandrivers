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


################################### DATA #######################################
# Read and inspect data

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