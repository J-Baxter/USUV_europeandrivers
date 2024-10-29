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
options(scipen = 6, digits = 4) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)
library(igraph)

# User functions
# Calculates pairwise hamming distance, then groups at a threshold difference
GroupSequences <- function(aln, snp_threshold = 0){
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
  
  if( any(hd_raw[lower.tri(hd_raw, diag = FALSE)] <= snp_threshold)){
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
  }else{
    print('all sequences above threshold.')
    groups <- tibble(tipnames = rownames(aln)) %>%
      rowid_to_column(., var = 'sequence_group')
  }
  
    
    out <- tibble(tipnames = rownames(aln)) %>%
      left_join(groups) %>%
      mutate(sequence_group = 
               ifelse(is.na(sequence_group), 
                      max(sequence_group, na.rm = T) + row_number() + 1, 
                      sequence_group))

  
  return(out)
}

############################################## DATA ################################################
nflg_alignment_files <- c('./2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_A.fasta',
                          './2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_B.fasta',
                          './2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_C.fasta',
                          './2024Oct20/alignments/subset_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_D.fasta')

nflg_subset_alignments <- lapply(nflg_alignment_files,
                                 read.dna,
                                 as.matrix = T,
                                 format = 'fasta')



# import metadata
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')

metadata_split <- lapply(nflg_subset_alignments, function(x) metadata %>% filter(tipnames %in% rownames(x)))

############################################## MAIN ################################################

# exclude Tempest outliers

to_exclude <- 'MK230893|Grivegnee/2017|turdus_merula|BE|2017-08'

nflg_subset_alignments %<>% 
  lapply(., function(x) x[!rownames(x) %in% to_exclude,])


# Subsampling 1 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across each cluster alignment

groups <- lapply(nflg_subset_alignments,
                  GroupSequences,
                 snp_threshold = 1) %>%
  
  setNames(., LETTERS[1:4]) %>%
  bind_rows(., .id = 'CLUSTER_GROUP')



subsample_1 <- metadata_split %>%
  
  # set names
  setNames(., LETTERS[1:4]) %>%
  
  bind_rows(., .id = 'CLUSTER_GROUP') %>%
  
  # join identity groups
  left_join(groups) %>%
  
  # establish groupings
  group_by(CLUSTER_GROUP,
           date_y,
           sequence_group) %>%
  
  slice_sample(n = 1) %>%
  
  ungroup() %>%
  group_split(CLUSTER_GROUP)


aln_subsampled <- mapply(function(aln, meta) aln[rownames(aln) %in% meta$tipnames,],
                         nflg_subset_alignments,
                         subsample_1)



############################################## WRITE ###############################################
filenames <- gsub('.fasta', '_subsampled.fasta', nflg_alignment_files)

mapply(write.dna,
       aln_subsampled,
       filenames,
       format = 'fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################
