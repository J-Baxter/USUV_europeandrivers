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
  hd_normalised <- dist.hamming(aln) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  groups <- which(hd_raw <= snp_threshold, 
                  arr.ind = TRUE) %>%
    dplyr::as_tibble(rownames = 'tipnames') %>% 
    filter(row !=col) %>%
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
nflg_alignment <- read.dna('./2024Dec02/alignments/USUV_2024Dec02_alldata_aligned_formatted_noFLI_NFLG.fasta',
                      as.matrix = T,
                      format = 'fasta')


# import metadata
metadata_nflg <- read_csv('./data/USUV_metadata_all_2024Dec02.csv') %>%
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


# Subsampling 2
# this uses the tip order from the MCC tree calculated using the previous subsample

# get list of tips from major european clusters on MCC tree
eu_clusters <- offspring(nflg_mcc, c(737, 624, 450, 852, 435), tiponly = TRUE) %>%
  lapply(., function(x) dplyr::select(as_tibble(nflg_mcc), c(label, node)) %>% 
           filter(node %in% x) %>%
           pull(label)) %>%
  set_names(c('A', 'B', 'C')) %>%
  lapply(., as_tibble) %>%
  bind_rows(,.id = 'cluster') %>%
  rename(tipnames = value)
  

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
          './2024Dec02/alignments/USUV_2024Dec02_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta')

write.nexus.data()

#write.dna(aln_subsampled_2,
       #   './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample2.fasta', 
        #  format = 'fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################
