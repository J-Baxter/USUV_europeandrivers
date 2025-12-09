################################################################################
## Script Name:        Calculate Bioclimatic Variables
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-11-24
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


################################### DATA #######################################
# Read and inspect data

################################### MAIN #######################################
# Main analysis or transformation steps

# BIO1 - Annual Mean Temperature
# Units: Degrees Celsius
# Data Inputs: The average temperature for each month

# BIO2 — Annual Mean Diurnal Range
# Units: Degrees Celsius
# Data Inputs: Monthly maximum and monthly minimum temperatures

# BIO3 - Isothermality
# Units: Percent
# Data Inputs: Results from equation 2 and equation 8


# BIO4 - Temperature Seasonality
# Units: Temperature (degrees Celsius)
# Data Inputs: The average temperature for each month

# BIO5 - 
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################