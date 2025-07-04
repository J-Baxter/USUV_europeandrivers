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

################################### DATA #######################################
# Read and inspect data
nflg_hipstr <- read.beast('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_constant_hipstr.tree')

nflg_ml <- read.newick('./2025Jun24/global_analysis/phylogenetic_placement/USUV_2025Jun24_alldata_aligned_formatted_noFLI_NFLG_subsampled.fasta.contree') %>% 
  phytools::midpoint_root()

subsample_clusterings <- read_csv('./2025Jun24/europe_clusters/subsample_clusterings.csv')

subsample_nomenclature <- read_delim('./2025Jun24/nomenclature/subsample/USUV_2025Jun24_labels.tsv') %>%
  rename(label = sample)

################################### MAIN #######################################
# Main analysis or transformation steps
mrca_nodes <- nflg_ml %>%

  # ML tree to tree-data object and add cluster/lineage assignments
  as_tibble() %>%
  left_join(subsample_clusterings %>% dplyr::select(label, dist_30)) %>%
  
  # Exclude problematic cluster I
  filter(dist_30 != 'I') %>%

  # Exclude all NAs and split by cluster
  drop_na(dist_30) %>%
  group_split(dist_30) %>%
  setNames(LETTERS[1:length(.)]) %>%
  
  # Find common ancestor nodes of clusters/lineages 
  lapply(., function(x) MRCA(nflg_ml, .node1 = x$label)) %>%
  enframe(name = 'cluster', value = 'root') %>%
  mutate(root = as.integer(root))


# visual checks
#nflg_ml %>%
  #left_join(mrca_nodes %>% rename(node = root)) %>%
  #ggtree() + 
  #geom_nodepoint(aes(colour = cluster))


# Find all descendants of common ancestor nodes
all_seq_clusters <- mrca_nodes %>%
  pull(root) %>%
  lapply(. , function(x) ml_tbl %>% 
                     filter(node %in% offspring(nflg_ml, x, type = 'tips')) %>%
                     pull(label)) %>%
  lapply(., as_tibble) %>%
  setNames(LETTERS[1:length(.)]) %>%
  bind_rows(, .id = 'cluster') %>%
  rename(label = value)

# visual checks
nflg_ml %>%
  left_join(all_seq_clusters) %>% 
  ggtree(.) + 
  geom_tippoint(aes(colour = cluster))
  

# Cluster X lineage comparison
subsample_clusterings %>% dplyr::select(label, dist_30) %>%
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
  left_join(all_seq_clusters) %>% 
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

# all nflg clusters
write_csv(all_seq_clusters, './2025Jun24/europe_clusters/all_nflg_cluster.csv')

#################################### END #######################################
################################################################################
