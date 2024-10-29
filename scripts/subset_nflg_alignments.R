####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2024-10-29
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(tidyverse)
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)

# User functions


############################################## DATA ################################################
nflg_alignment <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta',
                           as.matrix = T,
                           format = 'fasta')

metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')

clusters <- read_csv('./data/clusters_2024Oct29.csv')

metadata %<>%
  left_join(clusters)

############################################## MAIN ################################################

split_metadata <- metadata %>%
  
  # all europe NFLG allocated a cluster - exclude the rest for now
  drop_na(cluster) %>%
  
  group_split(cluster)

split_alignments <- lapply(split_metadata, function(x) nflg_alignment[rownames(nflg_alignment) %in% x$tipnames,])

############################################## WRITE ###############################################
filenames <- paste0('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG_',
                   c('A', 'B', 'C', 'D', 'E'),
                   '.fasta')

mapply(write.dna,
       split_alignments,
       filenames,
       format = 'fasta')

############################################## END #################################################
####################################################################################################
####################################################################################################