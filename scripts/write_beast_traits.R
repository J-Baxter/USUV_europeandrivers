################################################################################
## Script Name:        Write Trait files and Edit XMLs
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-07-11
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
library(ape)


################################### DATA #######################################
# Read and inspect data
metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')

nflg_cluster_subsampled_files <- list.files('./2025Jun24/alignments',
                                            pattern = 'NFLG_[A-Z]+_subsampled',
                                            full.names = T)

nflg_cluster_subsampled <- lapply(nflg_cluster_subsampled_files,
                                  read.dna,
                                  format = 'fasta',
                                  as.matrix = T)


partial_cluster_subsampled_files <- list.files('./2025Jun24/alignments',
                                               pattern = 'partial_[A-Z]+_subsampled',
                                               full.names = T)

partial_cluster_subsampled <- lapply(partial_cluster_subsampled_files,
                                     read.dna,
                                     format = 'fasta',
                                     as.matrix = T)

################################### MAIN #######################################
# Main analysis or transformation steps
metadata_for_beast <- metadata_with_concat %>%
  dplyr::select(tipnames, ends_with('id'), geocode_coords, location_precision) %>%
  mutate(geocode_coords = gsub('c\\(|\\,|\\)', '', geocode_coords)) %>%
  separate_wider_delim(geocode_coords, delim = ' ', names = c('long', 'lat')) %>%
  mutate(long = case_when(nuts3_id == 'ITH35' & location_precision == 'nuts3' ~ '12.55874', .default = long),
         lat = case_when(nuts3_id == 'ITH35' & location_precision == 'nuts3' ~ '45.60318', .default = lat)) %>%
  dplyr::select(-location_precision) %>%
  relocate(tipnames, nuts0_id, nuts1_id, nuts2_id, nuts3_id, lat, long) 


metadata_for_beast_nflg <- nflg_cluster_subsampled %>%
  lapply(., rownames) %>%
  lapply(., function(x) metadata_for_beast %>% filter( tipnames %in% x))

metadata_for_beast_partial <- partial_cluster_subsampled %>%
  lapply(., rownames) %>%
  lapply(., function(x) metadata_for_beast %>% filter( tipnames %in% x))


################################### OUTPUT #####################################
# Save output files, plots, or results
nflg_filename <- gsub('noFLI_alldata_aligned_formatted', 
                      '',
                      nflg_cluster_subsampled_files) %>%
  gsub('.fasta$', '_traits.txt', .) %>%
  gsub('alignments', 'europe_clusters', .) 

partial_filename <- gsub('noFLI_alldata_aligned_formatted',
                         '', 
                         partial_cluster_subsampled_files) %>%
  gsub('.fasta$', '_traits.txt', .) %>%
  gsub('alignments', 'europe_clusters', .) 

mapply(write_delim,
       metadata_for_beast_nflg,
       nflg_filename,
       delim = '\t',
       quote= 'needed')

mapply(write_delim,
       metadata_for_beast_partial,
       partial_filename,
       delim = '\t',
       quote= 'needed')

#################################### END #######################################
################################################################################