################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-09-10
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
# library(tidyverse)
library(magrittr)
library(ape)


################################### DATA #######################################
# Read and inspect data
unclocklike <- read_csv('./2025Jun24/europe_clusters/cluster_phylo_ml_2/partial_clusters_unclocklike_2.csv')


partial_cluster_alignment_master_files <- list.files(path = './2025Jun24/alignments',
                                                     pattern = 'USUV_2025Jun24_alldata_aligned_formatted_noFLI_partial_[:A-Z:]{1,4}_subsampled.fasta',
                                                     full.names = TRUE)

partial_cluster_alignment_working_files <- list.files('./2025Jun24/europe_clusters',
                                                      pattern = 'partial_[:A-Z:]{1,4}_subsampled.fasta$',
                                                      recursive = T,
                                                      full.names = T) 


################################### MAIN #######################################
# Main analysis or transformation steps
partial_cluster_alignments <- lapply(c(partial_cluster_alignment_master_files,
                                       partial_cluster_alignment_working_files),
                                     read.dna,
                                     as.matrix = T,
                                     format = 'fasta') 

partial_cluster_alignments_filtered <- lapply(partial_cluster_alignments, 
                                              function(x) x[!rownames(x) %in% unclocklike$label,])



################################### OUTPUT #####################################
# Save output files, plots, or results
filtered_filenames <- gsub('\\.fasta$', '_filtered\\.fasta',
                           )


mapply(write.FASTA,
       partial_cluster_alignments_filtered,
       c(partial_cluster_alignment_master_files,
         partial_cluster_alignment_working_files))  

#################################### END #######################################
################################################################################. 