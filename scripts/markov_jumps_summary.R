test <- read_delim('./2025Jun24/europe_clusters/workspace/nuts0_id.jumpHistory.log', delim = '\t', skip = 3) %>%
  head() %>%
  rowwise() %>%
  mutate(jumps = str_split(completeHistory_1, '\\}\\,\\{')) %>%
  unnest_longer(jumps) %>%
  mutate(jumps = str_remove(jumps, '\\{\\{|\\}\\}')) %>%
  separate_wider_delim(jumps, delim = ',', names = c('site', 'time', 'origin', 'destination'))