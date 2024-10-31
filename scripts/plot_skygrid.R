skygrid_b <- read_delim('./2024Oct20/alignments/subset_alignments/run/B_skygrid.txt') %>%
  separate_wider_delim(5, '\t', names = c('mean', 'median', 'upper', 'lower')) %>% 
  .[-1,] %>%
  rename(time = ...1) %>%
  mutate(across(c(1,5,6,7,8), .fns = ~as.numeric(.x)))

ggplot(skygrid_b) +
  geom_line(aes(x = time, y = median), colour = 'brown') + 
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = 'brown', alpha = 0.1)+
  scale_y_log10()+
  theme_minimal()