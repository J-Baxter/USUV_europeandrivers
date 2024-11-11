# Dispersal velocity drivers

velocity_resistance <- read_csv('./2024Oct20/alignments/concatenated_alignments/leastcostmodelresults_velocity.csv')
velocity_conductance <- read_csv('./2024Oct20/alignments/concatenated_alignments/leastcostmodelresults_velocity_conductance.csv')

# testing if the corresponding environmental variable as a factor repulsing or attracting tree nodes, respectively.

direction_files <- list.files(pattern = '(E|R)_Bayes') 

direction <-  list.files(pattern = '(E|R)_Bayes')   %>%
  lapply(., read.table, header = T) %>%
  setNames( list.files(pattern = '(E|R)_Bayes') ) %>%
  bind_rows(., .id = 'filename' ) %>%
  drop_na(BF) %>%
  as_tibble() %>%
  mutate(cluster = str_split_i(filename, '_', 2),
        var= gsub('USUV_[[:upper:]]_|_least-cost.*', '', filename) %>%
          #gsub('[[:punct:]]+', ' ', .) %>%
          str_to_lower(),
        statistic = gsub('.*direction_|_Bayes.*', '', filename)) %>%
  select(-1) %>%
  pivot_wider(names_from = statistic, values_from = BF) %>%
  rename(remain_in = E,
         leave_from = R)


#E = the mean of the environmental values extracted at the nodes’ position, 
#R = the proportion of branches for which the environmental value recorded at the oldest node
# position is higher than the environmental value recorded at the youngest node position. 

#While E measures the tendency of tree nodes to remain located in lower/higher environmental values, 
# R rather measures the tendency of lineages to disperse towards lower/higher environmental values.


direction_strong <- direction %>% filter(E>20 | R >20)