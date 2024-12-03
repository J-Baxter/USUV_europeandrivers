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

# Only non-published sequences
original_format %<>%
  filter(is.na(GenBank))

erasmus_data <- original_format %>% 
  dplyr::select(-c( 1,
                   DateInfo,
                   English,
                   family,
                   order,
                   Remark,
                   Surveillance,
                   CT)) %>%
  mutate(Country = case_when(!is.na(location) ~ 'Netherlands')) %>%
  
  unite(., coords, Lat, Long, sep = ', ') %>%
  mutate(Isolate = coalesce(UniqueID, ID),
         location = coalesce(location, coords)) %>%
  unite(., 'Geo_Location', postal.code, location, province, sep = ' ', na.rm = T) %>%
  dplyr::select(-c(ID, UniqueID, coords)) %>%
  rename(Collection_Date = DeathOrSampleDate,
         Host = Latin,
         Isolation_Source = Sample,
         Accession = GenBank) %>%
  
  # where multiple sequences from the sample bird exist (due to different isolation source), down-
  # sample. Preference for Throat_Cloacal_Swab samples
  group_by(Isolate) %>%
  filter(!(Isolation_Source != 'Throat_Cloacal_Swab' & n() > 1)) %>%
  ungroup()

erasmus_labels <- erasmus_data %>%
  unite(., label, Isolate, Collection_Date, sep = '|') %>%
  pull(label)

erasmus_sequences <- erasmus_data %>%
  pull(Sequence) %>%
  lapply(., function(x) str_split(x, '') %>% unlist() %>% as.matrix() %>% t()) %>%
  setNames(erasmus_labels) %>%
  as.DNAbin()



############################################## WRITE ###############################################

write_csv(erasmus_data %>% dplyr::select(-Sequence),
          './data/erasmus_data.csv')

write.FASTA(erasmus_sequences,
            './erasmus_data/erasmus_sequences.fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################