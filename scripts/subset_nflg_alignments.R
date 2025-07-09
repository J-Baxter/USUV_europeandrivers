################################################################################
## Script Name:        Split the NFLG alignment according to a csv of lineage/cluster designations+
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


################################### DATA #######################################
# Read and inspect data
nflg_alignment <- read.dna('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG.fasta',
                           as.matrix = T,
                           format = 'fasta')

metadata <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')
clusters <- read_csv('./2025Jun24/europe_clusters/all_nflg_cluster.csv')

metadata %<>%
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

split_alignments <- lapply(split_metadata_by_cluster, function(x) 
  nflg_alignment[rownames(nflg_alignment) %in% x$tipnames,])

################################### OUTPUT #####################################
# Save output files, plots, or results
filenames <- paste0('./2025Jun24/alignments/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_',
                    as.roman(1:length(split_alignments)),
                    '.fasta')

mapply(write.FASTA,
       split_alignments,
       filenames)
#################################### END #######################################
################################################################################