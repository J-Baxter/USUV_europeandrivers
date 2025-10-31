################################################################################
## Script Name:        Patristic Distance Figures
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-10-31
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
library(treeio)
library(TreeTools)
library(ape)
library(ade4)
library(adephylo)
library(ggtree)
library(ggsci)


InferClusters <- function(phylo,
                          n_threshold = 3, 
                          dist_threshold =50,
                          filter = TRUE,
                          metadata){
  
  stopifnot(class(phylo) == 'phylo')
  
  
  # Get subtrees from main tree and filter by minimum tip threshold
  subtrees <- subtrees(phylo) %>%
    .[lapply(., Ntip) %>% unlist() >= n_threshold]
  
  if(isTRUE(filter)){
    int <- lapply(subtrees, function(tree) metadata %>% filter(tipnames %in% tree$tip.label)) %>%
      lapply(., function(x) x %>% dplyr::pull(is_europe)) %>% 
      lapply(., function(x) all(x == 1)) %>%
      unlist() %>%
      which()
    
    subtrees <-  subtrees[int]
  }
  
  
  # Calculate patristic distance matrix for all subtrees
  subtree_patristic <- lapply(subtrees, adephylo::distTips) %>%
    lapply(., as.matrix) 
  
  # determine maximum patristic distance for each subtree
  subtree_maxpat <- lapply(subtree_patristic, max) %>% 
    unlist()
  
  # filter subtrees (upper bound patristic distance)
  subtree_below_threshold <-  subtrees[which(subtree_maxpat < dist_threshold)] 
  
  # Determine which subtrees are fully subsumed by other subtrees
  tips_per_subtree <- lapply(subtree_below_threshold, function(x) x$tip.label) 
  
  part_of <- vector(mode = "list", length = length(tips_per_subtree))
  
  for (i in 1:length(tips_per_subtree)){
    
    tips_per_subtree[[1]]
    
    for ( j in 1:length(tips_per_subtree)){
      
      if(all(tips_per_subtree[[i]] %in% tips_per_subtree[[j]]) & i!=j){
        
        part_of[[i]] <- c(part_of[[i]], j)
        
      }
    }
  }
  
  # Elements left null are the 'top tier' clusters 
  int_parent_nodes <- lapply(part_of, is.null) %>%
    unlist() %>% 
    which()
  
  clusters <- subtree_below_threshold[int_parent_nodes]
  
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  
  # Format output tibble containing tip names and cluster designation
  cluster_tips <- lapply(clusters, TipLabels) %>%
    setNames(LETTERS702[1:length(.)]) %>%
    lapply(., cbind.data.frame) %>%
    lapply(., as_tibble) %>%
    bind_rows(., .id = 'cluster') %>%
    rename(label = `X[[i]]`)
  
  return(cluster_tips)
  
}


################################### DATA #######################################
# Read and inspect data
nflg_hipstr <- read.beast('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_constant_hipstr.tree')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')%>%
  filter(tipnames %in% nflg_hipstr@phylo$tip.label) 


################################### MAIN #######################################
# Main analysis or transformation steps
cluster_verylong <- lapply(seq(5, 80, by = 5),
                           InferClusters, 
                           phylo = nflg_hipstr@phylo,
                           n_threshold = 2, 
                           filter = TRUE,
                           metadata_in_tree) %>%
  setNames(as.character(seq(5, 80, by = 5))) %>%
  bind_rows(., .id = 'distance_threshold') 

cluster_wide <- cluster_verylong %>%
  pivot_wider(values_from = cluster,
              names_from = distance_threshold,
              names_prefix = 'dist_')


# Calculate within/between Europe clade tips
europe_tips <- metadata_in_tree %>% 
  filter(is_europe == T) %>%
  pull(tipnames)

tmp <- dist_matrix %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(cluster_wide) 

tmp_2 <- tmp %>%
  group_by(pair_id) %>%
  dplyr::select(pair_id ,distance, name , label,
                dist_5,dist_15, dist_25, dist_35, dist_45, dist_55, dist_65, dist_75) %>%
  mutate(across(
    starts_with("dist_"),
    ~ case_when(
      n_distinct(.x) == 1 ~ 'Within Europe: Within Clade',
      TRUE ~ 'Within Europe: Between Clade'
    ),
    .names = "type_{str_extract(.col, '\\\\d+')}"
  ))



# Elbow curve plotting the number of clusters agains the patristic distance
# threshold
cluster_verylong %>%
  summarise(n = n_distinct(cluster), .by = distance_threshold) %>%
  ggplot() + 
  geom_vline(aes(xintercept = 35), linetype = 'dashed', colour = 'red')+
  geom_point(aes(y = n, x = as.numeric(distance_threshold))) +
  scale_y_continuous('Clusters (n)') + 
  scale_x_continuous('Patristic Distance Threshold', expand = c(0,0)) + 
  my_theme 


# 2. Tree and patristic distance distributions for original/
most_recent_date <- '2024-10-12'

nflg_hipstr %>% 
  left_join(cluster_wide) %>%
  ggtree(mrsd = most_recent_date) + 
  scale_y_reverse() +
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_35)) + 
  scale_color_brewer(palette = 'Dark2', direction = -1)



# 3. patristic distance historgrams at different thresholds
tmp_2 %>%
  pivot_longer(c(starts_with('type') ),
               names_to = 'cluster_threshold',
               values_to = 'cluster_type') %>%
  mutate(cluster_threshold = gsub('type_', '', cluster_threshold) %>%
           as.numeric()) %>%
  select(pair_id, distance, label, starts_with('cluster')) %>%
  ggplot() +
  geom_histogram(aes(x = distance, 
                     fill = cluster_type, 
                     colour = cluster_type),
                 alpha = 0.5,
                 binwidth = 2.5,
                 position="identity" )+
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL) + 
  scale_fill_d3(palette = 'category20',
                alpha = 0.99,
                name= NULL) +
  scale_x_continuous('Patristic Distance',
                     expand = c(0,0))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic()+
  theme(legend.position = 'bottom') +
  facet_wrap(~ cluster_threshold, 
             ncol = 1,
             strip.position = 'right',
             scales = 'free_x')

################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################