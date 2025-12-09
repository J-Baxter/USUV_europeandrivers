################################################################################
## Script Name:       Calculate Wavefront, Diffusion Coefficients and IBD
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
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
library(tidyverse)
library(magrittr)
library(seraphim)
library(giscoR)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(diagram)
library(lubridate)
library(treeio)
ExtractTrees <- function(tree,
                         localTreesDirectory,
                         mostRecentSamplingDatum){
  tab = postTreeExtractions(post_tre=tree, mostRecentSamplingDatum)
  
  i_in_file = list.files(localTreesDirectory) %>%
    str_extract_all(., '\\d+') %>%
    as.numeric() %>%
    max(na.rm = T)
  
  i <- ifelse(i_in_file >= 1, i_in_file + 1, 1)
  
  write.csv(tab,
            paste0(localTreesDirectory, "/TreeExtractions_", i, ".csv"), 
            row.names=F, 
            quote=F)
}


CalcSpreadStatistics <- function(treefile, 
                         localTreesDirectory,
                         mostRecentSamplingDatum,
                         burnIn = 0, 
                         randomSampling = FALSE, 
                         nberOfTreesToSample = 100,
                         coordinateAttributeNam = "location"){
  
  # Extracting spatio-temporal information embedded in trees with the 
  # "postTreeExtractions" function (advised)
  trees = readAnnotatedNexus(treefile)
  dir.create(localTreesDirectory, showWarnings=F)
  
  lapply(trees,
         ExtractTrees, 
         localTreesDirectory = localTreesDirectory,
         mostRecentSamplingDatum = mostRecentSamplingDatum)


  # 2. Estimation of several dispersal statistics
  nberOfExtractionFiles = 100
  timeSlices = 100
  onlyTipBranches = FALSE
  showingPlots = FALSE
  outputName = gsub("^(.*?/[^/_]+)_.*$", "\\1", treefile) # first _ after last /
  nberOfCores = 8
  slidingWindow = 1
  
  spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, 
                   showingPlots, outputName, nberOfCores, slidingWindow)
  
}

################################### DATA #######################################

################################### MAIN #######################################
# Main analysis or transformation steps
CalcSpreadStatistics('./2025Jun24/europe_clusters/Partial_III/USUV_2025Jun24_partial_III_SRD06_fixed_SG_joint_500.trees',
                     localTreesDirectory = './2025Jun24/europe_clusters/Partial_III/Extracted_trees',
                     mostRecentSamplingDatum = 2024.669)

CalcSpreadStatistics('./2025Jun24/europe_clusters/Partial_V/USUV_2025Jun24_partial_V_SRD06_fixed_SG_joint_500.trees',
                     localTreesDirectory = './2025Jun24/europe_clusters/Partial_V/Extracted_trees',
                     mostRecentSamplingDatum = 2023.564)

CalcSpreadStatistics('./2025Jun24/europe_clusters/Partial_VI/USUV_2025Jun24_partial_VI_SRD06_fixed_SG_joint_500.trees',
                     localTreesDirectory = './2025Jun24/europe_clusters/Partial_VI/Extracted_trees',
                     mostRecentSamplingDatum = 2024.861)

CalcSpreadStatistics('./2025Jun24/europe_clusters/Partial_VII/USUV_2025Jun24_partial_VII_SRD06_fixed_SG_joint_500.trees',
                     localTreesDirectory = './2025Jun24/europe_clusters/Partial_VII/Extracted_trees',
                     mostRecentSamplingDatum = 2016.658)

CalcSpreadStatistics('./2025Jun24/europe_clusters/Partial_VIII/USUV_2025Jun24_partial_VIII_SRD06_fixed_SG_joint_500.trees',
                     localTreesDirectory = './2025Jun24/europe_clusters/Partial_VIII/Extracted_trees',
                     mostRecentSamplingDatum = 2019.000)

################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################