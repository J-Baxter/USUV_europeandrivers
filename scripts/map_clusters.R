################################################################################
## Script Name:        Map clusters and lineages to all NFLG data
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
library(magrittr)
library(treeio)
library(TreeTools)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(phangorn)


MapClusters <- function(new_tree, cluster_tbl){
  
  tree_tbl <- as_tibble(new_tree)
  
  mrca_nodes <- new_tree %>%
    
    # ML tree to tree-data object and add cluster/lineage assignments
    as_tibble() %>%
    left_join(cluster_tbl %>% dplyr::select(label, dist_35)) %>%
    
    # Exclude problematic cluster I
    #filter(dist_35 != 'I') %>%
    
    # Exclude all NAs and split by cluster
    drop_na(dist_35) %>%
    group_split(dist_35) %>%
    setNames(LETTERS[1:length(.)]) %>%
    
    # Find common ancestor nodes of clusters/lineages 
    lapply(., function(x) MRCA(new_tree, .node1 = x$label)) %>%
    enframe(name = 'cluster', value = 'root') %>%
    mutate(root = as.integer(root))
  
  
  # Find all descendants of common ancestor nodes
  all_seq_clusters <- mrca_nodes %>%  
    pull(root) %>%
    lapply(. , function(x) tree_tbl %>% 
             filter(node %in% offspring(new_tree, x, type = 'tips')) %>%
             pull(label)) %>%
    lapply(., as_tibble) %>%
    setNames(LETTERS[1:length(.)]) %>%
    bind_rows(, .id = 'cluster') %>%
    rename(label = value)
  
  return(all_seq_clusters)
}


################################### DATA #######################################
# Read and inspect data

# Previously inferred data
# subsampled beast tree
nflg_hipstr <- read.beast('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_constant_hipstr.tree')
subsample_clusterings <- read_csv('./2025Jun24/europe_clusters/subsample_clusterings.csv')
subsample_nomenclature <- read_delim('./2025Jun24/nomenclature/subsample/USUV_2025Jun24_labels.tsv') %>%
  rename(label = sample)

# New data
# all NFLG sequences (ML)
nflg_ml <- read.newick('./2025Jun24/europe_clusters/nflg_phyloplacement/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG.fasta.contree') %>% 
  phytools::midpoint_root()

# subsampled NFLG and all concat/partial sequences (ML)
partial_ml <- read.newick('./2025Jun24/europe_clusters/partial_seq_phyloplacement/USUV_2025Jun24_alldata_aligned_formatted_subsampled_incpartial.fasta.treefile') %>% 
  phytools::midpoint_root()


################################### MAIN #######################################
# Main analysis or transformation steps
all_nflg_clusters <- MapClusters(nflg_ml, subsample_clusterings)
all_partial_clusters <- MapClusters(partial_ml, subsample_clusterings) %>%
  left_join(metadata_with_concat %>% dplyr::select(label = tipnames, generegion_nflg))

all_clusters <-  bind_rows(all_nflg_clusters,
                           all_partial_clusters) %>%
  slice_sample(n =1, by = label)

# visual checks
partial_ml %>%
  left_join(all_partial_clusters) %>% 
  mutate(generegion_nflg = as.character(generegion_nflg)) %>%
  ggtree(.) + 
  geom_tippoint(aes(colour = cluster, shape = generegion_nflg))


# Cluster X lineage comparison
subsample_clusterings %>% dplyr::select(label, dist_35) %>%
  left_join(subsample_nomenclature) %>%
  dplyr::select(-1) %>%
  distinct() %>%
  view()

# visualise (ugly)
lineage_colours <- c('A' = '#a4243b',
                     'B' = '#be776b',
                     'C' = '#d8c99b',
                     'D' = '#d8b06c', 
                     # 'E' =  '#d8973c', 
                     'E' = '#cb7d36', 
                     'F' = '#bd632f', 
                     'G' = '#273e47')

nflg_hipstr %>%
  left_join(all_clusters) %>% 
  ggtree(mrsd = '2024-10-12') +
  geom_cladelab(data = level_0_tbl %>% rename(lineage_label = lineage) %>% 
                  mutate(level_0 = gsub('\\..*', '', lineage_label)),
                mapping = aes(node = root_node, 
                              label = lineage_label, 
                              colour = level_0),
                offset = -32,
                offset.text = -17,
                fontsize = 3, 
                #textcolour = NA,
                #size = 0.01,
                align = T) + 
  scale_colour_manual(values = lineage_colours, guide = 'none')+
  new_scale_colour() +
  
  # Label level 1
  geom_cladelab(data = level_1_tbl %>% rename(lineage_label = lineage) %>% 
                  mutate(level_0 = gsub('\\..*', '', lineage_label)),
                mapping = aes(node = root_node, 
                              label = lineage_label, 
                              colour = level_0),
                offset = -20,
                offset.text = -17,
                # size = 0.01,
                align = TRUE, 
                fontsize = 3, 
                #textcolour = NA,
                #barsize = 15,
  ) + 
  scale_colour_manual(values = lineage_colours, guide = 'none') +
  new_scale_colour() + 
  geom_tippoint(aes(colour = cluster))


################################### OUTPUT #####################################
# Save output files, plots, or results
# all clusters
write_csv(all_clusters, './2025Jun24/europe_clusters/all_clusters.csv')
#################################### END #######################################
################################################################################