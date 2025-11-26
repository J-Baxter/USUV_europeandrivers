library(purrr)


CalulcateLogLikelihood <- function(dataframe){
  dataframe %<>%
    mutate(type = str_to_lower(type) %>% 
             gsub('clade$|within europe: | ', '', .)) 

  summary_stats <- dataframe %>%
    summarise(mu = mean(distance), 
              sigma = sd(distance),
              n_clust = n_distinct(pick(6)) ,.by = type) %>%
    pivot_wider(names_from = type, 
                values_from = c(mu, sigma))
  
  dataframe %<>% 
    bind_cols(summary_stats)
  
  
  log_likelihood <- dataframe %>%
    split(~type) %>% 
    
    map_at("between",  
           ~ mutate(.x, loglikelihood = dnorm(distance,
                                              mean = mu_between, 
                                              sd = sigma_between, 
                                              log = T))) %>%
    map_at('within',  
           ~ mutate(.x, loglikelihood = dnorm(distance, 
                                              mean = mu_within, 
                                              sd = sigma_within, 
                                              log = T))) %>%
    bind_rows() %>%
    pull(loglikelihood) %>%
    sum()
  
  return(log_likelihood)
}

CalulcateLogLikelihood(pdist30_clades)
CalulcateLogLikelihood(pdist40_clades)
CalulcateLogLikelihood(pdist50_clades)
CalulcateLogLikelihood(pdist60_clades)



teststat <- -2 * (CalulcateLogLikelihood(pdist60_clades) -CalulcateLogLikelihood(pdist50_clades))
pchisq(teststat, df = 1, lower.tail = FALSE)
