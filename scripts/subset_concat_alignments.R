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

clusters <- read_csv('./data/clusters_2024Oct29.csv')

concat_metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv') %>%
  left_join(clusters)

ml_concat <- read.tree('./2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_withconcatenated.fasta.contree')

concat_metadata %<>%
  left_join(clusters)



############################################## MAIN ################################################

# first root ml tree by south african sequences
rooted_concat <- root(ml_concat, outgroup =  "EU074021|NA|NA|ZA|1958" , resolve.root = TRUE)


# Determine the MRCA of clusters, as determined by NFLG sequences
cluster_mrcas <- rooted_concat %>%
  as.treedata() %>%
  left_join(concat_metadata,
            by = join_by(label == tipnames)) %>%
  as_tibble() %>%
  drop_na(cluster) %>%
  group_split(cluster) %>%
  lapply(., pull, 'node') %>%
  lapply(., getMRCA, phy = rooted_concat)


cluster_offspring <- lapply(cluster_mrcas, offspring, .data =  rooted_concat, tiponly = TRUE) %>%
  lapply(., function(x) rooted_concat$tip.label[x]) %>% 
  lapply(., as_tibble) %>%
  setNames(LETTERS[1:5]) %>%
  bind_rows(., .id = 'cluster') %>%
  rename(tipnames = value)


# Subset alignments and update metadata
cluster_metadata <- concat_metadata %>%
  rows_patch(., cluster_offspring, by='tipnames') 

cluster_metadata_split <- cluster_metadata %>%
  group_split(cluster)

cluster_concat_alignments <- lapply(cluster_metadata_split, function(x) concat_alignment[rownames(concat_alignment) %in% x$tipnames,])


############################################## WRITE ###############################################
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
  geom_hilight(node=1683, fill="#1B9E77") + 
  geom_hilight(node=1739, fill= "#D95F02" )+
  geom_hilight(node=2031, fill="#7570B3")+
  geom_hilight(node=2323, fill= "#E7298A")+
  geom_hilight(node=2321, fill= "#66A61E")
  


filenames <- paste0('./2024Oct20/alignments/concatenated_alignments//USUV_2024Oct20_alldata_aligned_formatted_noFLI_concat_',
                    c('A', 'B', 'C', 'D', 'E'),
                    '.fasta')

mapply(write.dna,
       cluster_concat_alignments[-6],
       filenames,
       format = 'fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################