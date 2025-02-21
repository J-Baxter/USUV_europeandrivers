Accordingly, we could then 
# determine log-likelihood of observing the patristic distance matrix for ‘k’ clusters at patristic 
# distance threshold ‘t’:  We selected the optimal threshold according to log-likelihood ratio tests
# at successive patristic distance thresholds. 

list(pdist30_clades,
     pdist40_clades,
     pdist50_clades,
     pdist60_clades) %>%
  
  # At each patristic distance threshold, we calculated the mean and standard deviation of the 
  # within-lineage and between-lineage patristic distance distributions. 
  lapply(., function(x) x %>%
           summarise(mean = mean(distance), sd = sd(distance), n_clust = n_distinct(pick(6)) ,.by = type)) %>%
  setNames(c(30,40,50,60)) %>%
  bind_rows(, .id = 'threshold') 
  
  # determine the log-likelihood of observing the patristic distance matrix for ‘k’ clusters at patristic 
  # distance threshold ‘t’: 

  
  
distance_matrix <- as.matrix(patristic_distances)


pair_loglikelihood <- pmap()
n <- nrow(distance_matrix)
log_likelihood <- 0

# Iterate over all pairs of sequences
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    distance <- distance_matrix[i, j]
    
    if (cluster_assignments[i] == cluster_assignments[j]) {
      # Same cluster: Use within-cluster parameters
      likelihood_ij <- normal_pdf(distance, mu_within, sigma_within)
    } else {
      # Different clusters: Use between-cluster parameters
      likelihood_ij <- normal_pdf(distance, mu_between, sigma_between)
    }
    
    # Add the log-likelihood of this pair
    log_likelihood <- log_likelihood + log(likelihood_ij)
  }
}



# Calculate the log-likelihood of a patristic distance matrix under a clustering hypothesis
calculate_log_likelihood <- function(distance_matrix, cluster_assignments, 
                                     mu_within, sigma_within, 
                                     mu_between, sigma_between) {
  n <- nrow(distance_matrix)
  log_likelihood <- 0
  
  # Iterate over all pairs of sequences
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      distance <- distance_matrix[i, j]
      
      if (cluster_assignments[i] == cluster_assignments[j]) {
        # Same cluster: Use within-cluster parameters
        likelihood_ij <- normal_pdf(distance, mu_within, sigma_within)
      } else {
        # Different clusters: Use between-cluster parameters
        likelihood_ij <- normal_pdf(distance, mu_between, sigma_between)
      }
      
      # Add the log-likelihood of this pair
      log_likelihood <- log_likelihood + log(likelihood_ij)
    }
  }
  
  return(log_likelihood)
}

# Example usage
# Example patristic distance matrix (symmetric matrix)
distance_matrix <- matrix(c(0, 2, 5, 6,
                            2, 0, 4, 5,
                            5, 4, 0, 3,
                            6, 5, 3, 0), 
                          nrow = 4, byrow = TRUE)

# Example cluster assignments
cluster_assignments <- c(1, 1, 2, 2)

# Parameters for within-cluster and between-cluster distances
mu_within <- 2
sigma_within <- 0.5
mu_between <- 5
sigma_between <- 1

# Calculate the log-likelihood
log_likelihood <- calculate_log_likelihood(distance_matrix, cluster_assignments, 
                                           mu_within, sigma_within, 
                                           mu_between, sigma_between)
print(log_likelihood)
