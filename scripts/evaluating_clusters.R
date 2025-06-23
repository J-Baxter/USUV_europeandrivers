################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-06-02
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


################################### DATA #######################################
# Read and inspect data
nflg_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/lineage_level_0_taxa/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_lineage0taxa_mcc.tree')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025May22.csv')%>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

################################### MAIN #######################################
# Main analysis or transformation steps

cluster_verylong <- lapply(seq(5, 80, by = 5),
                           InferClusters, 
                           phylo = nflg_mcc@phylo,
                           n_threshold = 2, 
                           filter = TRUE,
                           metadata_in_tree) %>%
  setNames(as.character(seq(5, 80, by = 5))) %>%
  bind_rows(., .id = 'distance_threshold') 


# Elbow method
cluster_verylong %>%
  summarise(n = n_distinct(cluster), .by = distance_threshold) %>%
  ggplot() + 
  geom_point(aes(y = n, x = as.numeric(distance_threshold)))


cluster_wide <- cluster_verylong %>%
  pivot_wider(values_from = cluster,
              names_from = distance_threshold,
              names_prefix = 'dist_')

most_recent_date <- '2024-10-12'

nflg_mcc %>% 
  left_join(cluster_wide) %>%
  ggtree(mrsd = most_recent_date) + 
  scale_y_reverse() +
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_30))

# Example data
# dist_matrix: symmetric matrix of pairwise distances
# clusters: vector of cluster assignments

# Convert dist object to matrix if needed
dist_matrix <- adephylo::distTips(nflg_mcc@phylo) %>% as.matrix()

CalcLikelihood <- function(dist_mat, clusters, ignore.na = TRUE){
  n <- nrow(dist_mat)
  
  # Create all pairs (i < j)
  pairs <- which(upper.tri(dist_mat), arr.ind = TRUE)
  
  # Assign labels to each pair
  same_cluster <- clusters[pairs[,1]] == clusters[pairs[,2]]
  
  if (ignore.na == TRUE){
    valid_pairs <- !is.na(clusters[pairs[,1]]) & !is.na(clusters[pairs[,2]])
    
    # Filter pairs and recalculate
    pairs_valid <- pairs[valid_pairs, , drop = FALSE]
    same_cluster <- clusters[pairs_valid[,1]] == clusters[pairs_valid[,2]]
    
    within_dists <- dist_matrix[pairs_valid][same_cluster]
    between_dists <- dist_matrix[pairs_valid][!same_cluster]
  }
  
  # Estimate Gaussian parameters
  mu_within <- mean(within_dists)
  sigma_within <- sd(within_dists)
  mu_between <- mean(between_dists)
  sigma_between <- sd(between_dists)
  
  # Compute log-likelihood
  log_likelihood <- function(d, same) {
    if (same) {
      dnorm(d, mean = mu_within, sd = sigma_within, log = TRUE)
    } else {
      dnorm(d, mean = mu_between, sd = sigma_between, log = TRUE)
    }
  }
  
  # Apply to all pairs
  log_liks <- mapply(log_likelihood,
                     d = dist_matrix[pairs_valid],
                     same = same_cluster)
  
  total_log_likelihood <- sum(log_liks)
  
  return(total_log_likelihood)
}

logLik_bm <- ace(cluster_list[[1]], nflg_mcc@phylo, type = "discrete")$loglik

cluster_list <-tibble(label = rownames(dist_matrix)) %>%
  left_join(cluster_wide) %>%
  pivot_longer(-1, values_to = 'cluster', names_to = 'threshold') %>%
  select(-label) %>%
  group_split(threshold) %>%
  lapply(., function(x) x %>% pull(cluster))


log_likelihoods <- lapply(cluster_list, function(x) CalcLikelihood(dist_matrix, x)) %>%
  flatten_dbl()

k <- cluster_list %>% lapply(., function(x) length(unique(x))) %>%
  flatten_dbl()

bic <- -2 * log_likelihoods + k * log(n)
aic <- 2 * k - 2 * log_likelihoods

tibble(bic = bic,
       aic = aic,
          threshold = seq(5, 80, by = 5)) %>%
  ggplot() + 
  geom_point(aes(x = threshold, y = aic))
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################








pdist_clades <- lapply()
  
  
  test <- function(x, dist_mat){
    
  }
  
  
tmp <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
    
  left_join(cluster_wide) 


tmp_2 <- tmp %>%
  group_by(pair_id) %>%
  mutate(across(
    starts_with("dist_"),
    ~ case_when(
      n_distinct(.x) == 1 ~ 'Within Europe: Within Clade',
      TRUE ~ 'Within Europe: Between Clade'
    ),
    .names = "type_{str_extract(.col, '\\\\d+')}"
  ))


tmp_3 <- tmp_2 %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, starts_with('type')), names_from = 'name', values_from =  c(label, starts_with('dist_'))) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b))))
  


dists <- tmp_3 %>%
  dplyr::select(c(starts_with('type'), distance)) %>%
  mutate(across(starts_with('type'), .fns = ~  gsub('within europe\\: ', '', str_to_lower(.x)))) %>%
  pivot_longer(cols = starts_with('type'), names_to = 'threshold', values_to = 'group') %>%
  drop_na(group, threshold) %>%
  group_by(threshold, group) %>%
  reframe(stats = glance(fitdistr(distance, 'normal')))
  

tmp_data <- tmp_3 %>%
  dplyr::select(c(starts_with('type'), distance)) %>%
  mutate(across(starts_with('type'), .fns = ~  gsub('within europe\\: ', '', str_to_lower(.x)))) %>%
  pivot_longer(cols = starts_with('type'), names_to = 'threshold', values_to = 'group') %>%
  filter(threshold == 'type_35') %>%
  filter(group == 'between clade') %>%
  pull(distance)

fit <- fitdistr(tmp_data, 'normal')
glance(fit)
ggplot(tmp_3) +
  geom_histogram(aes(x = distance, 
                     #y = after_stat(ndensity),
                     fill = type_35, colour = type_35),
                 alpha = 0.5,
                 binwidth = 5,
                 position="identity" )+
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL) + 
  scale_fill_d3(palette = 'category20',
                alpha = 0.99,
                name= NULL) +
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic()+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))


