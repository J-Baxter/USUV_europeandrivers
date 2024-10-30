####################################################################################################
####################################################################################################
## Script name: edit_duplicated_isolates
##
## Purpose of script: to check for duplicate isolates in the metadata and either a) downsample if 
## at least one full genome is present or b) concatenate disjoint sequences into one.
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
metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')

alignment <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI.fasta',
                      format = 'fasta',
                      as.matrix = T)

############################################## MAIN ################################################
duplicate_isolates <- metadata %>%
  drop_na(sequence_isolate) %>%
  group_by(sequence_isolate, date_ymd) %>%
  filter(!any(generegion_nflg == 1)) %>%
  mutate(is_duplicate_isolate = n()) %>%
  ungroup() %>%
  filter(is_duplicate_isolate > 1) %>%
  group_split(sequence_isolate) 


duplicate_seqs <- duplicate_isolates %>%
  lapply(., function(x) alignment[rownames(alignment) %in% x$tipnames,]) 


concat_seqnames_tbl <- lapply(duplicate_seqs, rownames) %>%
  lapply(., as_tibble)%>%
  bind_rows() %>%
  rename(old_seqnames = value) %>%
  mutate(tipnames  = gsub('^[^|]*\\|', 'NA\\|', old_seqnames))


concat_seqnames <- concat_seqnames_tbl %>%
  pull(tipnames) %>%
  unique()

concat_alignments <- duplicate_seqs %>%
  lapply(., as.list) %>%
  lapply(.,as.character) %>%
  lapply(., function(x) lapply(x, function(y) gsub('-', NA, y))) %>%
  lapply(., function(x) coalesce(!!!x)) %>%
  lapply(., function(x) replace_na(x, '-')) %>%
  do.call(rbind, .) %>%
  set_rownames(concat_seqnames) %>%
  as.DNAbin()


# Create new metadata rows for concatenated sequences and down sample repeat NFLG sequences
updated_metadata <- metadata %>%
  drop_na(sequence_isolate) %>%
  group_by(sequence_isolate) %>%
  mutate(is_duplicate_isolate = n()) %>%
  ungroup() %>%
  filter(is_duplicate_isolate > 1) %>%
  group_by(sequence_isolate) %>%
  mutate(across(-any_of(c("generegion_NS5_9000_9600",
                          "generegion_NS5_10042_10312",
                          "eneregion_env_1003_1491",
                          "generegion_nflg")), .fns = ~coalesce(.x))) %>%
  mutate(old_seqnames = tipnames) %>%
  rows_update(concat_seqnames_tbl) %>%
  slice_sample(n=1) %>%
  ungroup() %>%
  dplyr::select(-old_seqnames)
   

excluded_tipnames <-  metadata %>% filter(tipnames %in% c(updated_metadata$tipnames, concat_seqnames_tbl$old_seqnames)) %>%
  pull(tipnames) 


new_metadata <- metadata %>%
  filter(!tipnames %in% excluded_tipnames) %>%
  bind_rows(updated_metadata)

new_alignment <- alignment %>%
  rbind(concat_alignments) %>%
  .[rownames(.) %in% new_metadata$tipnames,]




############################################## WRITE ###############################################

write.dna(new_alignment,
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_withconcatenated.fasta',
          format = 'fasta')


write_csv(new_metadata,
          './data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv')


############################################## END #################################################
####################################################################################################
####################################################################################################