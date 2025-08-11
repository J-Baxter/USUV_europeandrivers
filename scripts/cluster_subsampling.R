################################################################################
## Script Name:        Downsample alignments for within-Europe analysis (main)
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
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


SubSampleClusterAlignments <- function(aln_list, data, exclude_seqs = NULL){
  
  cluster_tbl <- lapply(aln_list, rownames) %>%
    set_names(1:length(.)) %>%
    lapply(., as_tibble) %>%
    bind_rows(., .id = 'cluster') %>%
    mutate(cluster = as.numeric(cluster)) %>%
    rename(tipnames = value)
  
  hd_groups <- lapply(aln_list,
                      GroupSequences,
                      snp_threshold = 20) %>%
    
    set_names(1:length(.)) %>%
    bind_rows(., .id = 'cluster')%>%
    mutate(cluster = as.numeric(cluster)) 
  
  if(!missing(exclude_seqs)){
    old_nrow <- nrow(data)
    data <- data %>% filter(!tipnames %in% exclude_seqs)
  }
  
  # Details stratified subsampling aiming to 'thin' the tree by selecting one 
  # sequence per NUTS3, per year-month. 
  # This is applied across each cluster alignment (clusters < 5 sequences excluded)
  subsampled <- data %>%
    mutate(date_quarter = quarter(as.Date(paste0(date_ym, '-01')))) %>%
    left_join(cluster_tbl) %>%
    arrange(cluster) %>%
    left_join(hd_groups, by = c('cluster', 'tipnames')) %>%
    drop_na(cluster) %>%
    group_by(cluster) %>%
    
    # Only analyse clusters with five or more sequences
    #filter(n() >= 5) %>%
    
    group_by(date_ym,
             #date_quarter,
             eurostat_polygon,
             sequence_group,
             .add = TRUE) %>%
    
    slice_sample(n = 1) %>%
    
    ungroup() %>%
    group_split(cluster)
  
  test <- subsampled %>%
    bind_rows() %>%
    filter(tipnames %in% exclude_seqs)
  
  stopifnot(nrow(test) == 0)
  # summarise(n_seqs = n(), n_eurostat = n_distinct(eurostat_polygon), .by = cluster)
  
  subsampled_clusters <- lapply(subsampled, function(x) x %>% 
                                  pull(cluster) %>% 
                                  unique()) %>%
    unlist() %>%
    as.numeric()
  
  out <- mapply(
    function(aln, x) aln[rownames(aln) %in% x$tipnames,],
    aln_list[names(aln_list) %in% subsampled_clusters],
    subsampled)
  
  return(out)
  
}
################################### DATA #######################################
# Read and inspect data
nflg_cluster_alignment_files <- list.files(path = './2025Jun24/alignments',
                                           pattern = 'USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_[[A-Z]]*',
                                           full.names = TRUE)

partial_cluster_alignment_files <- list.files(path = './2025Jun24/alignments',
                                           pattern = 'USUV_2025Jun24_alldata_aligned_formatted_noFLI_partial_[[A-Z]]*',
                                           full.names = TRUE)

nflg_cluster_alignments <- lapply(nflg_cluster_alignment_files,
                                  read.dna,
                                  as.matrix = T,
                                  format = 'fasta') %>%
  set_names(1:length(.))

partial_cluster_alignments <- lapply(partial_cluster_alignment_files,
                                  read.dna,
                                  as.matrix = T,
                                  format = 'fasta') %>%
  set_names(1:length(.))

metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')

tempest_check <- read_csv('./2025Jun24/europe_clusters/cluster_phylo_ml/tempest_check.csv')
################################### MAIN #######################################
# Main analysis or transformation steps
# 1. Exclude unclocklike sequences
to_exclude <- tempest_check %>%
  drop_na(notes) %>%
  filter(grepl('\\|', notes)) %>%
  pull(notes) %>%
  # separate out multiple entries
  str_split(., ', ') %>%
  unlist()


# 2. Subsample
nflg_cluster_subsampled <- SubSampleClusterAlignments(nflg_cluster_alignments,
                                                      metadata_with_concat,
                                                      exclude_seqs = to_exclude)

partial_cluster_subsampled <- SubSampleClusterAlignments(partial_cluster_alignments, 
                                                         metadata_with_concat,
                                                         exclude_seqs = to_exclude)

# 3. Keep only alignments that are viable (i.e partial_sequences > 10)
viable <-which(unlist(lapply(partial_cluster_subsampled, nrow) >= 10))


################################### OUTPUT #####################################
# Save output files, plots, or results
filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_',
                    as.roman(names(nflg_cluster_subsampled[viable])),
                    '_subsampled.fasta')

mapply(write.FASTA,
       nflg_cluster_subsampled[viable],
       filenames)  

partial_cluster_filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_partial_',
                                    as.roman(names(partial_cluster_subsampled[viable])),
                                    '_subsampled.fasta')

mapply(write.FASTA,
       partial_cluster_subsampled[viable],
       partial_cluster_filenames)  

  
#################################### END #######################################
################################################################################