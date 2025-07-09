################################################################################
## Script Name:        Downsample NFLGs for within-Europe analysis
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


################################### DATA #######################################
# Read and inspect data
nflg_cluster_alignment_files <- list.files(path = './2025Jun24/alignments',
                                           pattern = 'USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_[[A-Z]]*',
                                           full.names = TRUE)

nflg_cluster_alignments <- lapply(nflg_cluster_alignment_files,
                                  read.dna,
                                  as.matrix = T,
                                  format = 'fasta') %>%
  set_names(1:length(.))

metadata <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')


################################### MAIN #######################################
# Main analysis or transformation steps
clusters <- lapply(nflg_cluster_alignments, rownames) %>%
  set_names(1:length(.)) %>%
  lapply(., as_tibble) %>%
  bind_rows(., .id = 'cluster') %>%
  mutate(cluster = as.numeric(cluster)) %>%
  rename(tipnames = value)

hd_groups <- lapply(nflg_cluster_alignments,
                 GroupSequences,
                 snp_threshold = 20) %>%
  
  set_names(1:length(.)) %>%
  bind_rows(., .id = 'cluster')%>%
  mutate(cluster = as.numeric(cluster)) 
  

# Details stratified subsampling aiming to 'thin' the tree by selecting one 
# sequence per NUTS3, per year-month. 
# This is applied across each cluster alignment (clusters < 5 sequences excluded)
subsampled <- metadata %>%
  mutate(date_quarter = quarter(as.Date(paste0(date_ym, '-01')))) %>%
  left_join(clusters) %>%
  arrange(cluster) %>%
  left_join(hd_groups, by = c('cluster', 'tipnames')) %>%
  drop_na(cluster) %>%
  group_by(cluster) %>%
  
  # Only analyse clusters with five or more sequences
  filter(n() >= 5) %>%
  
  group_by(date_ym,
           #date_quarter,
           eurostat_polygon,
           sequence_group,
           .add = TRUE) %>%
  
  slice_sample(n = 1) %>%
  
  ungroup() %>%
  group_split(cluster)

  summarise(n_seqs = n(), n_eurostat = n_distinct(eurostat_polygon), .by = cluster)

subsampled_clusters <- lapply(subsampled, function(x) x %>% 
                                pull(cluster) %>% 
                                unique()) %>%
  unlist() %>%
  as.numeric()

nflg_cluster_subsampled <- mapply(
  function(aln, x) aln[rownames(aln) %in% x$tipnames,],
  nflg_cluster_alignments[names(nflg_cluster_alignments) %in% subsampled_clusters],
  subsampled)

subsampled %>%
  summarise(n_seqs = n(), .by = c(sequence_group, eurostat_polygon, cluster)) %>%
  ggplot(aes(x = as.character(sequence_group), y = as.character(eurostat_polygon), fill = n_seqs)) +
  geom_raster() +
  facet_wrap(.~cluster, scales = 'free')

metadata %>%
  left_join(clusters) %>%
  left_join(hd_groups, by = c('cluster', 'tipnames')) %>%
  drop_na(cluster) %>%
  summarise(n_seqs = n(), n_eurostat = n_distinct(eurostat_polygon), .by = cluster)


################################### OUTPUT #####################################
# Save output files, plots, or results
filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_',
                    as.roman(names(nflg_cluster_subsampled)),
                    '_subsampled.fasta')

mapply(write.FASTA,
       nflg_cluster_subsampled,
       filenames)  

  
#################################### END #######################################
################################################################################


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



############################################## DATA ################################################
############################################## MAIN ################################################

# exclude Tempest outliers

to_exclude <- 'MK230893|Grivegnee/2017|turdus_merula|BE|2017-08'

nflg_subset_alignments %<>% 
  lapply(., function(x) x[!rownames(x) %in% to_exclude,])

############################################## END #################################################
####################################################################################################
####################################################################################################
