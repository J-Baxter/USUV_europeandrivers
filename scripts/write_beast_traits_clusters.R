####################################################################################################
####################################################################################################
## Script name: 
##
## Purpose of script: create beast metadata files for within-europe lineages. these include discrete 
## trait phylogeography at NUTS0 and NUTS1 levels; also continuous phylogeography (with option for 
## heterogenous prior on uncertain regions)
##
## Date created: 2024-10-30
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)

# User functions


############################################## DATA ################################################
aln_subsampled_files <- c("./2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_A_subsampled.fasta",
                          "./2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_B_subsampled.fasta",
                          "./2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_C_subsampled.fasta",
                          "./2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_D_subsampled.fasta") 

nflg_subsampled_alignments <- lapply(aln_subsampled_files,
                                     read.dna,
                                     as.matrix = T,
                                     format = 'fasta')

metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')
metadata_split <- lapply(nflg_subsampled_alignments, function(x) metadata %>% filter(tipnames %in% rownames(x)))


############################################## MAIN ################################################
beast_data <- metadata_split %>% 
  
  # set names
  setNames(., LETTERS[1:4]) %>%
  
  bind_rows(., .id = 'CLUSTER_GROUP') %>%
  
  # spread geocode coords to lat lon
  mutate(geocode_coords = gsub('c|\\(|\\)|,', '', geocode_coords)) %>%
  separate_wider_delim(geocode_coords, ' ', names =  c('lon', 'lat')) %>%
  mutate(across(c(lat, lon), .fns = ~ as.numeric(.x))) %>%

  
  # Select columns for BEAST
  dplyr::select(
    tipnames,
    lat,
    lon,
    nuts0_id,
    nuts1_id,
    CLUSTER_GROUP) %>%
  group_split(CLUSTER_GROUP, .keep = FALSE)



############################################## WRITE ###############################################
beast_filename <- gsub('_alldata_aligned_formatted', '', aln_subsampled_files) %>%
  gsub('.fasta$', '_traits.txt', .)

mapply(write_delim,
       beast_data,
       beast_filename,
       delim = '\t',
       quote= 'needed')





############################################## END #################################################
####################################################################################################
####################################################################################################