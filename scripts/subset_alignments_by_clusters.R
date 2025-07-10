################################################################################
## Script Name:        Split alignments according to a csv of cluster designations+
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-07-04
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
library(tidyverse)
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)

SplitAlignmentByCluster <- function(alignment, data){
  split_metadata_by_cluster <- data %>%
    drop_na(cluster) %>%
    group_split(cluster)
  
  # Sanity check - each sequence included only once (and all sequences are present)
  stopifnot(lapply(split_metadata_by_cluster, nrow) %>% unlist() %>% sum() == 
              nrow(data %>% drop_na(cluster)))
  
  out <- lapply(split_metadata_by_cluster, function(x) 
    alignment[rownames(alignment) %in% x$tipnames,])
  
  return(out)
  
}

################################### DATA #######################################
# Read and inspect data

# all NFLG sequences
nflg_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG.fasta',
                           as.matrix = T,
                           format = 'fasta')

partial_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_incpartial.fasta',
                                      as.matrix = T,
                                      format = 'fasta')

metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')
clusters <- read_csv('./2025Jun24/europe_clusters/all_clusters.csv')

metadata_with_concat %<>%
  left_join(clusters,
            by = join_by('tipnames' == 'label'))


################################### MAIN #######################################
# Main analysis or transformation steps
split_metadata_by_cluster <- metadata %>%
  
  # all europe NFLG allocated a cluster - exclude the rest for now
  drop_na(cluster) %>%
  group_split(cluster)

# Sanity check - each sequence included only once (and all sequences are present)
stopifnot(lapply(split_metadata_by_cluster, nrow) %>% unlist() %>% sum() == 
            nrow( metadata %>% drop_na(cluster)))

split_nflg_alignments <- SplitAlignmentByCluster(nflg_alignment, 
                                                 metadata_with_concat)

split_partial_alignments <- SplitAlignmentByCluster(partial_alignment, 
                                                    metadata_with_concat)

################################### OUTPUT #####################################
# Save output files, plots, or results
nflg_filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_',
                    as.roman(1:length(split_nflg_alignments)),
                    '.fasta')

mapply(write.FASTA,
       split_nflg_alignments,
       nflg_filenames)


partial_filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_partial_',
                         as.roman(1:length(split_partial_alignments)),
                         '.fasta')

mapply(write.FASTA,
       split_partial_alignments,
       partial_filenames)
#################################### END #######################################
################################################################################