####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2024-10-28
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions


############################################## DATA ################################################
nflg_mcc<- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_mcc.tree')
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')
patristic_distance_clusters <- read_csv('./2024Oct20/nflg_patristic_distance_clusters.csv')


metadata_in_tree <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv') %>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(nflg_mcc@phylo)) #sanity check

############################################## MAIN ################################################

beast_data <- metadata_in_tree %>% 
  # spread geocode coords to lat lon
  mutate(geocode_coords = gsub('c|\\(|\\)|,', '', geocode_coords)) %>%
  separate_wider_delim(geocode_coords, ' ', names =  c('lon', 'lat')) %>%
  mutate(across(c(lat, lon), .fns = ~ as.numeric(.x))) %>%
  
  # format binary trait
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                               .default = 'not_europe')) %>%
  
  left_join(patristic_distance_clusters,
            by = join_by(tipnames == label)) %>%
  
  # Select columns for BEAST
  dplyr::select(
    tipnames,
    lat,
    lon,
    is_europe,
    starts_with('dist_')) %>%
  
  mutate(across(starts_with('dist'), .fns = ~ case_when(grepl('KU760915', tipnames) ~ 'X',
                                                       grepl('not_europe', is_europe) ~ 'source',
                                                       .default = .x)))


############################################## WRITE ###############################################
beast_filename <-paste('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_traits.txt')

write_delim(beast_data,
            beast_filename,
            delim = '\t',
            quote= 'needed')

  


############################################## END #################################################
####################################################################################################
####################################################################################################