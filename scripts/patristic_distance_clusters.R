####################################################################################################
####################################################################################################
## Script name: Cluster inference 
##
## Purpose of script: Estimate clusters based on the maximum patristic distance of sub trees to 
## qualitatively compare import scenarios
##
## Date created: 2024-10-29
##
##
########################################## SYSTEM OPTIONS ##########################################
options(ignore.negative.edge=TRUE)


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(treeio)
library(TreeTools)
library(ape)
library(ade4)
library(adephylo)
library(ggtree)

# User functions
InferClusters <- function(phylo, n_threshold = 3, dist_threshold =50, filter = TRUE, metadata){
  
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



############################################## DATA ################################################
nflg_ca <- read.beast('./2025Feb10/global_analysis/beast_plain/USUV_2025Feb10_NFLG_SRD06_relaxLn_constant_mcc.tree')
#nflg_mcc <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_mcc.tree')

metadata <- read_csv('./data/USUV_metadata_all_2025Feb10.csv')


############################################## MAIN ################################################
InferClusters(phylo = nflg_ca@phylo,
              metadata, 
              n_threshold = 3, 
              filter = TRUE,
              dist_threshold = 50)

# Get clusters with maximum patristic distances of 30, 40, 50 and 60
cluster_long <- lapply(c(30, 40, 50, 60),
       InferClusters, 
       phylo = nflg_ca@phylo,
       n_threshold = 2, 
       filter = TRUE,
       metadata) %>%
  setNames(c('30', '40', '50', '60')) %>%
  bind_rows(., .id = 'distance_threshold') 


cluster_wide <- cluster_long %>%
  pivot_wider(values_from = cluster,
              names_from = distance_threshold,
              names_prefix = 'dist_')

# Test plot
most_recent_date <- metadata %>%
  filter(tipnames %in% nflg_ca@phylo$tip.label) %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym
  max(na.rm = TRUE)

nflg_ca %>% 
  left_join(cluster_wide) %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_30)) # Facet not possible - must run as independent plots then 
# combined in cowplot (idea: one column of tree with corresponding distribution opposite)

############################################## WRITE ###############################################
write_csv(cluster_wide,
          './2024Oct20/nflg_patristic_distance_clusters.csv')

############################################## END #################################################
####################################################################################################
####################################################################################################