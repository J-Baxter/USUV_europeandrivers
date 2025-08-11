################################################################################
## Script Name:        Create Alignment of Partial & Concatenated Sequences
## Purpose:            Identify partial sequences, create concatenated alignment
#                      if multiple sequences for the same isolate exist, and 
#                      strip FLI alignments to ENV and NS5 CDS.
## Author:             James Baxter
## Date Created:       2025-07-09
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
full_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted.fasta',
                           as.matrix = T,
                           format = 'fasta')

nflg_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG.fasta',
                                 as.matrix = T,
                                 format = 'fasta')

nflg_subsampled_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta',
                                      as.matrix = T,
                                      format = 'fasta')

metadata <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')

################################### MAIN #######################################
# Main analysis or transformation steps
#### Originally partial sequences #### 
partial_tipnames <-  metadata %>% 
  filter(generegion_nflg == 0) %>%
  pull(tipnames)

partial_seq <- full_alignment[rownames(full_alignment) %in% partial_tipnames,]


#### Create gene-only sequences from FLI (ENV (976-2475) & NS5 (7684-10398)) ###
fli_tipnames <- metadata %>%
  filter(drop_fli == 1) %>%
  pull(tipnames)

fli_alignment <- full_alignment[rownames(full_alignment) %in% fli_tipnames, ] 
cols_to_keep <- c((976-97):(2475-97), (7684-97):(10398-79))
cols_to_mask <- setdiff(seq_len(ncol(fli_alignment)), cols_to_keep)
fli_alignment[, cols_to_mask] <- as.raw(240)


#### Concatenated alignments from multiple partial seqs #### 
duplicate_isolates <- metadata %>%
  drop_na(sequence_isolate) %>%
  group_by(sequence_isolate, date_ymd) %>%
  filter(!any(generegion_nflg == 1)) %>%
  filter(n() > 1) %>%
  ungroup()

# Split alignment into different isolates
duplicate_seqs <- duplicate_isolates %>%
  group_split(sequence_isolate) %>%
  lapply(., function(x) full_alignment[rownames(full_alignment) %in% x$tipnames,]) 

# Edit tipnames
concat_seqnames_tbl <- lapply(duplicate_seqs, rownames) %>%
  lapply(., as_tibble)%>%
  bind_rows() %>%
  rename(old_seqnames = value) %>%
  mutate(tipnames  = gsub('^[^|]*\\|', 'NA\\|', old_seqnames))


concat_seqnames <- concat_seqnames_tbl %>%
  pull(tipnames) %>%
  unique()

# Merge concatenated alignments
concat_alignment <- duplicate_seqs %>%
  lapply(., as.list) %>%
  lapply(.,as.character) %>%
  lapply(., function(x) lapply(x, function(y) gsub('-', NA, y))) %>%
  lapply(., function(x) coalesce(!!!x)) %>%
  lapply(., function(x) replace_na(x, '-')) %>%
  do.call(rbind, .) %>%
  set_rownames(concat_seqnames) %>%
  as.DNAbin()


#### join all alignments #### 
joint_align <- rbind(
  nflg_alignment,
  fli_alignment,
  concat_alignment,
  partial_seq[!rownames(partial_seq) %in% concat_seqnames_tbl$old_seqnames,]) 
  
joint_align <- joint_align[unique(rownames(joint_align)),]

## subsampled alignment for phylogenetic placement
joint_subsampled_align <- rbind(
  nflg_subsampled_alignment,
  fli_alignment,
  concat_alignment,
  partial_seq[!rownames(partial_seq) %in% concat_seqnames_tbl$old_seqnames,]) 

joint_subsampled_align <- joint_subsampled_align[unique(rownames(joint_subsampled_align)),]


#### Update metadata #### 
duplicate_metadata <- duplicate_isolates %>%
  group_by(sequence_isolate) %>%
  mutate(across(-any_of(c("generegion_NS5_9000_9600",
                          "generegion_NS5_10042_10312",
                          "eneregion_env_1003_1491",
                          "generegion_nflg")), .fns = ~coalesce(.x))) %>%
  mutate(old_seqnames = tipnames,
         sequence_accession = paste(sequence_accession, collapse = ', ')) %>%
  rows_update(concat_seqnames_tbl) %>%
  slice_sample(n=1) %>%
  ungroup() %>%
  dplyr::select(-old_seqnames)

excluded_tipnames <-  metadata %>% 
  filter(tipnames %in% c(duplicate_metadata$tipnames, concat_seqnames_tbl$old_seqnames)) %>%
  pull(tipnames) 

new_metadata <- metadata %>%
  filter(!tipnames %in% excluded_tipnames) %>%
  bind_rows(duplicate_metadata)


################################### OUTPUT #####################################
# Save output files, plots, or results

write.FASTA(joint_subsampled_align,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_subsampled_incpartial.fasta')


write.FASTA(joint_align,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_incpartial.fasta')

write_csv(new_metadata,
          './data/USUV_metadata_2025Jun24_withconcatenated.csv')
#################################### END #######################################
################################################################################