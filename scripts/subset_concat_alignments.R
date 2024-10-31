####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2024-10-31
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
concat_alignment <- read.dna('./2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_withconcatenated.fasta',
                           as.matrix = T,
                           format = 'fasta')

ml_concat <- read.tree('./2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_withconcatenated.fasta.contree')
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')

clusters <- read_csv('./data/clusters_2024Oct29.csv')

metadata %<>%
  left_join(clusters)



############################################## MAIN ################################################

# first root ml tree by south african sequences
rooted_concat <- root(ml_concat, outgroup =  "EU074021|NA|NA|ZA|1958" , resolve.root = TRUE)


# Determine the MRCA of clusters, as determined by NFLG sequences
rooted_concat %>%
  as.treedata() %>%
  left_join(metadata,
            by = join_by(label == tipnames)) %>%
  as_tibble() %>%
  drop_na(cluster) %>%
  group_split(cluster) %>%
  lapply(., pull, 'node') %>%
  lapply(., getMRCA, phy = rooted_concat)



# Plot ML tree with labelled clades
rooted_concat %>%
  as.treedata() %>%
  left_join(metadata,
            by = join_by(label == tipnames)) %>%
  ggtree()+
  #geom_tippoint(aes(colour = cluster, shape = is_europe == 1))+
  scale_shape_manual(values = c("TRUE" = 18),
                     'New Sequences') +
  
  new_scale_colour()+
  geom_hilight(node=1683, fill="gold") + 
  geom_hilight(node=1739, fill="purple")+
  geom_hilight(node=2031, fill="green")+
  geom_hilight(node=2323, fill="blue")+
  geom_hilight(node=2321, fill="red")+
  

# Subset alignments and update metadata

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################