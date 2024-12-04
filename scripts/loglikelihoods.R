library(purrr)
summary_stats <- list(pdist30_clades,
     pdist40_clades,
     pdist50_clades,
     pdist60_clades) %>%
  
  # At each patristic distance threshold, we calculated the mean and standard deviation of the 
  # within-lineage and between-lineage patristic distance distributions. 
  lapply(., function(x) x %>%
           summarise(mean = mean(distance), sd = sd(distance), n_clust = n_distinct(pick(6)) ,.by = type)) %>%
  setNames(c(30,40,50,60)) %>%
  bind_rows(, .id = 'threshold') %>%
  mutate(type = gsub('Within Europe: | ', '', type) %>% str_to_lower()) %>%
  pivot_wider(names_from = type, values_from = c(mean, sd))


# Function to calculate the log-likelihood of a distance using cumulative probabilities
log_likelihood_distance <- function(dist, same_cluster, mu_within, sigma_within, mu_between, sigma_between) {
  if (same_cluster) {
    # Within-cluster distances: use within-cluster parameters
    density <- dnorm(dist, mean = mu_within, sd = sigma_within, log = T)
  } else {
    # Between-cluster distances: use between-cluster parameters
    density <- dnorm(dist, mean = mu_between, sd = sigma_between, log = T)
  }
  
  return()
}

# Calculate log-likelihood using pmap
calculate_log_likelihood_pmap <- function(distance_matrix, cluster_assignments, 
                                          mu_within, sigma_within, 
                                          mu_between, sigma_between) {
  # Get all unique pairs of indices from the distance matrix
    n <- nrow(distance_matrix)
    pairs <- expand.grid(i = 1:n, j = 1:n) %>% 
      filter(i < j)  # Keep only upper triangular pairs
  
  # Add distances and clustering information
  pairs <- pairs %>%
    mutate(distance = mapply(function(i, j) distance_matrix[i, j], i, j),
           same_cluster = cluster_assignments[i] == cluster_assignments[j]) %>%
    as_tibble() %>%
    select(distance, same_cluster)
  
  # Use pmap to calculate log-likelihoods for all pairs
  log_likelihoods <- apply(pairs, 1,
    log_likelihood_distance,  mu_within, sigma_within, mu_between, sigma_between, simplify =F) %>% 
    unlist()
)
  
  # Sum log-likelihoods
  total_log_likelihood <- sum(log_likelihoods)
  return(total_log_likelihood)
}

# Example usage
# Example patristic distance matrix (symmetric matrix)
distance_matrix <- as.matrix(patristic_distances)

distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cluster_wide$label, 
                                          colnames(distance_matrix) %in% cluster_wide$label]

# Example cluster assignments
cluster_assignments <- cluster_wide %>% select(c(label, dist_30)) %>%
  .[match(colnames(distance_matrix_subset), .$label),] %>%
  pull(dist_30)


length()
# Parameters for within-cluster and between-cluster distances
mu_within <- 15.0
sigma_within <- 5.76
mu_between <- 87.4
sigma_between <- 43.2

A tibble: 4 × 6
threshold n_clust mean_betweenclade mean_withinclade sd_betweenclade sd_withinclade
<chr>       <int>             <dbl>            <dbl>           <dbl>          <dbl>
  1 30             10              87.4             15.0            43.2           5.76
2 40              6              92.6             17.1            41.1           7.21
3 50              5              93.1             17.5            40.9           7.62
4 60              3             118.              28.3            26.3          15.0 


# Calculate the log-likelihood
log_likelihood <- calculate_log_likelihood_pmap(distance_matrix_subset, 
                                                cluster_assignments, 
                                                mu_within,
                                                sigma_within, 
                                                mu_between, 
                                                sigma_between)

