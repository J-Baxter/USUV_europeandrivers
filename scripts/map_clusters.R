####################################################################################################
####################################################################################################
## Script name: Map clusters to all NFLG sequences
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
library(treeio)
library(TreeTools)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(phangorn)


# User functions


############################################## DATA ################################################
nflg_traits <- read.beast('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits_mcc.tree')
nflg_mcc <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_mcc.tree')

nflg_ml <- read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta.treefile') %>% 
  phytools::midpoint_root()


metadata_in_beast_tree <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv') %>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

stopifnot(nrow(metadata_in_beast_tree ) == Ntip(nflg_mcc@phylo)) #sanity check

############################################## MAIN ################################################
ml_tbl <- as_tibble(nflg_ml)


root_nodes <- c('A' = 595,
                'B' = 751,
                'C' = 839,
                'D' = 578,
                'E' = 581)


# clusters inferred based on dist 50 threshold
clusters <- lapply(root_nodes, function(x) ml_tbl %>% 
                     filter(node %in% offspring(nflg_ml, x)) %>%
                     pull(label)) %>%
  lapply(., as_tibble) %>%
  bind_rows(, .id = 'cluster') %>%
  rename(tipnames = value)

############################################## WRITE ###############################################
write_csv(clusters,
      './data/clusters_2024Oct29.csv')
############################################## END #################################################
####################################################################################################



