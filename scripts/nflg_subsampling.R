####################################################################################################
####################################################################################################
## Script name: nflg downsampling
##
## Purpose of script:
##
## Date created: 2024-10-28
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)
library(phangorn)
library(igraph)


# User functions
# Calculates pairwise hamming distance, then groups at a threshold difference
GroupSequences <- function(aln, snp_threshold = 0){
  
  require(phangorn)
  require(igraph)
  
  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln) 
  }else{
    aln_formatted <- aln
  }
  
  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  groups <- which(hd_raw <= snp_threshold, 
                  arr.ind = TRUE) %>%
    dplyr::as_tibble(rownames = 'tipnames') %>% 
    filter(row != col) %>%
    dplyr::select(-tipnames) %>%
    
    # Infer network from HDs
    igraph::graph_from_data_frame(.,
                                  directed = F) %>%
    
    components() %>%
    getElement('membership') %>%
    stack() %>%
    as_tibble() %>%
    mutate(ind = as.numeric(as.character(ind))) %>%
    mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
    dplyr::select(c(tipnames, values)) %>%
    dplyr::distinct() %>%
    dplyr::rename(sequence_group = values)
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}

############################################## DATA ################################################
nflg_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG.fasta',
                      as.matrix = T,
                      format = 'fasta')


# import metadata
metadata_nflg <- read_csv('./data/USUV_metadata_all_2025Jun24.csv') %>%
  filter(generegion_nflg == 1)
  

############################################## MAIN ################################################

# exclude sequences that are too short or contain too many ambiguities
to_exclude <- metadata_nflg %>%
  filter(sequence_ambig > 0.05) %>%
  pull(tipnames) %>%
  
  # exclude Tempest outliers - need to redo
  c(.,'MK230893|Grivegnee/2017|turdus_merula|BE|2017-08'#,
    #'OQ414983|TV193/2019|culex_pipiens|RO|2019-07',
    #'KT445930|13-662|culex_modestus|CZ|2013-08-07'
    )
  
  
nflg_alignment %<>% .[!rownames(.) %in% to_exclude,]


# Subsampling 1 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across the entire 'global' dataset.

groups <- GroupSequences(nflg_alignment, 5) 

#groups %>%
  #summarise(n = n(), .by = sequence_group) %>%
 # pull(n) %>%
  #hist()


subsample_1 <- metadata %>%
  
  # include only sequences in alignment
  filter(tipnames %in% rownames(nflg_alignment)) %>%
  
  # join identity groups
  left_join(groups) %>%
  
  # establish groupings
  group_by(sequence_group,
           date_y,
           nuts0_id
           ) %>%
  
  slice_sample(n = 1)


aln_subsampled <- nflg_alignment[rownames(nflg_alignment) %in% subsample_1$tipnames,]


# Subsampling for sensitivity analysis
n_noneurope <- subsample_1 %>% 
  filter(is_europe == 0) %>%
  nrow()

n_europe <- subsample_1 %>% 
  filter(is_europe == 1) %>%
  nrow()

noneurope_subsamples <- subsample_1 %>% 
  filter(is_europe == 0) 

# 10% of downsampled europe sequences n = 45
groups_10 <- GroupSequences(nflg_alignment, 28) 
europe_10percent_subsample <- subsample_1 %>% 
  ungroup() %>%
  filter(is_europe == 1) %>%
  select(-sequence_group) %>%
  left_join(groups_10) %>%
  group_by(sequence_group,  nuts0_id) %>% # ngroups = 43
  slice_min(date_dec)
europe_10percent_subsample


# 25% of downsampled europe sequences n = 113
groups_25 <- GroupSequences(nflg_alignment, 21) 
europe_25percent_subsample <- subsample_1 %>% 
  ungroup() %>%
  filter(is_europe == 1) %>%
  select(-sequence_group) %>%
  left_join(groups_25) %>%
  group_by(sequence_group, date_y, nuts0_id) %>% # ngroups = 123
  slice_sample(n = 1)
europe_25percent_subsample


# 50% of downsampled europe sequences n = 227
groups_50 <- GroupSequences(nflg_alignment, 14) 
europe_50percent_subsample <- subsample_1 %>% 
  ungroup() %>%
  filter(is_europe == 1) %>%
  select(-sequence_group) %>%
  left_join(groups_50) %>%
  group_by(sequence_group, date_y, nuts1_id) %>% # ngroups = 230
  slice_sample(n = 1)
europe_50percent_subsample


# 75% of downsampled europe sequences n = 340
groups_75 <- GroupSequences(nflg_alignment, 8) 
europe_75percent_subsample <- subsample_1 %>% 
  ungroup() %>%
  filter(is_europe == 1) %>%
  select(-sequence_group) %>%
  left_join(groups_75) %>%
  group_by(sequence_group,  date_y , nuts0_id) %>% # ngroups = 357
  slice_sample(n = 1)



p75_subsample <- europe_75percent_subsample %>%
  bind_rows(noneurope_subsamples)

p75_subsample_aln <- nflg_alignment[rownames(nflg_alignment) %in% p75_subsample$tipnames,]


p50_subsample <- europe_50percent_subsample %>%
  bind_rows(noneurope_subsamples)

p50_subsample_aln <- nflg_alignment[rownames(nflg_alignment) %in% p50_subsample$tipnames,]


p25_subsample <- europe_25percent_subsample %>%
  bind_rows(noneurope_subsamples)

p25_subsample_aln <- nflg_alignment[rownames(nflg_alignment) %in% p25_subsample$tipnames,]


p10_subsample <- europe_10percent_subsample %>%
  bind_rows(noneurope_subsamples)

p10_subsample_aln <- nflg_alignment[rownames(nflg_alignment) %in% p10_subsample$tipnames,]

# this uses the tip order from the MCC tree calculated using the previous subsample

# get list of tips from major european clusters on MCC tree
#eu_clusters <- offspring(nflg_mcc, c(737, 624, 450, 852, 435), tiponly = TRUE) %>%
  #lapply(., function(x) dplyr::select(as_tibble(nflg_mcc), c(label, node)) %>% 
    #       filter(node %in% x) %>%
    #       pull(label)) %>%
  #set_names(c('A', 'B', 'C')) %>%
 # lapply(., as_tibble) %>%
  #bind_rows(,.id = 'cluster') %>%
  #rename(tipnames = value)
  

#subsample_2_noneu <- metadata %>%
  
  # include only sequences in alignment
  #filter(tipnames %in% rownames(nflg_alignment)) %>%
  #filter(is_europe == 0) %>%
  
 # left_join(groups) %>%
  
  # establish groupings
  #group_by(nuts1_id,
           #date_y) %>%
  #slice_sample(n = 1)

#subsample_2_eu <- metadata %>%
  
  # include only sequences in alignment
  #filter(tipnames %in% rownames(nflg_alignment)) %>%
  
  #left join clusters
  #left_join(eu_clusters) %>%
  
  # keep only main cluster 
 # filter(!is.na(cluster)) %>%
  
  #group_by(cluster,
          # is_europe) %>%

  # Sample first and last sequence from reassortant
 # slice(which(date_y == max(date_y, na.rm = T) | date_y == min(date_y, na.rm = T)), 
  #      .preserve = T) %>%
 # group_by(date_y, .add = TRUE) %>%
  #slice_sample(n = 1)
  
#subsample_2 <- bind_rows(subsample_2_noneu, subsample_2_eu)

#aln_subsampled_2 <- nflg_alignment[rownames(nflg_alignment) %in% subsample_2$tipnames,]


############################################## WRITE ###############################################
write.FASTA(aln_subsampled,
          './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta')

write.FASTA(p75_subsample_aln,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p75.fasta')

write.FASTA(p50_subsample_aln,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p50.fasta')

write.FASTA(p25_subsample_aln,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p25.fasta')

write.FASTA(p10_subsample_aln,
            './2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p10.fasta')

subsample_1 %>% 
  ungroup() %>%
  select(tipnames, is_europe) %>%
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                                   .default = 'not_europe')) %>%
  write_delim(.,
              './2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample1_traits.txt',
              delim = '\t',
              quote= 'needed')


p75_subsample %>% 
  ungroup() %>%
  select(tipnames, is_europe) %>%
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                               .default = 'not_europe')) %>%
  write_delim(.,
              './2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p75_traits.txt',
              delim = '\t',
              quote= 'needed')


p50_subsample %>% 
  ungroup() %>%
  select(tipnames, is_europe) %>%
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                               .default = 'not_europe')) %>%
  write_delim(.,
              './2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p50_traits.txt',
              delim = '\t',
              quote= 'needed')


p25_subsample %>% 
  ungroup() %>%
  select(tipnames, is_europe) %>%
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                               .default = 'not_europe')) %>%
  write_delim(.,
              './2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p25_traits.txt',
              delim = '\t',
              quote= 'needed')

p10_subsample %>% 
  ungroup() %>%
  select(tipnames, is_europe) %>%
  mutate(is_europe = case_when(is_europe == 1 ~'europe',
                               .default = 'not_europe')) %>%
  write_delim(.,
              './2025Jun24/global_analysis/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsample_p10_traits.txt',
              delim = '\t',
              quote= 'needed')


############################################## END #################################################
####################################################################################################
####################################################################################################
