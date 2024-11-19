####################################################################################################
####################################################################################################
## Script name: Extract dutch alignment and metadata from csv
##
## Purpose of script: as above
##
## Date created: 2024-11-18
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
original_format <- read_csv('./erasmus_data/240422_LiveDeadWildCapBirds_NL_16_22_USUVSeq.csv')


############################################## MAIN ################################################

# alignment with unique IDs
erasmus_labels <- original_format %>%
  unite(., label, ID, DeathOrSampleDate, sep = '|') %>%
  pull(label)

erasmus_sequences <- original_format %>%
  pull(Sequence) %>%
  lapply(., function(x) str_split(x, '') %>% unlist() %>% as.matrix() %>% t()) %>%
  setNames(erasmus_labels) %>%
  as.DNAbin()

erasmus_data <- original_format %>% 
  dplyr::select(-Sequence)

############################################## WRITE ###############################################

write_csv(erasmus_data,
          './erasmus_data/erasmus_data.csv')

write.FASTA(erasmus_sequences,
            './erasmus_data/erasmus_sequences.fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################