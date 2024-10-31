#playing with ggtree
library(tidyverse)
library(treeio)
library(ggtree)
library(ggmcmc)


GetNodeAges <- function(mrsd, nodeheight){
  # Convert the date string to a datetime object
  initial_date <-  as.Date(mrsd, format="%Y-%m-%d")
  
  # Calculate the number of days to subtract
  days_in_a_year <- 365.25  # Considering leap years
  days_to_subtract <- nodeheight * days_in_a_year
  
  # Subtract the days from the initial date
  final_date <- initial_date - days_to_subtract
  
  return(final_date)
  
}
all_europe <-  read.beast('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits_mcc.tree')
all_logs <- beastio::readLog('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits.log') %>%
  ggs()
tmrca <- all_logs %>% filter(Parameter == 'age.root.')


# treestat
treestat <- read_delim('~/Downloads/test.txt') %>%
  pivot_longer(., -1, names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = gsub("Height ", "", parameter) %>%
           as.numeric())

node_mrca <- treestat %>%
  filter(parameter == 116)

nflg_data <- data.frame(
  node = c(166, 197, 215, 141, 125, 118), 
  name = c("EU3", "EU2", 'EU1', 'AF3.3', 'AF3.2', 'AF3.1'))

# Define the shading intervals - age root for tmrca. Not sure yet on other nodes
shading_intervals <- seq(1850, 2040, by = 10)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2040)

# Create density plot
ggtree(all_europe, mrsd="2023-09-23") + 
  theme_tree2() +
  geom_tippoint(aes(color=collection_regionname))  +
  scale_colour_brewer(palette = 'Paired')+
  
  #Background
  annotate("rect", 
           xmin = shading_intervals[seq(1, length(shading_intervals), 2)], 
           xmax = shading_intervals_end[seq(1, length(shading_intervals_end), 2)], 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
  # Plot tmrca density
  geom_density(data = tmrca, 
               aes(x = value, y = (..count../Ntip(all_europe))),
               bounds = GetNodeAges("2023-09-23",  c(70.5184216919521, 179.175463927162)) %>% decimal_date() %>% sort(),
               fill = "blue", 
               alpha = 0.5, 
               inherit.aes = FALSE) +
  
  geom_vline(xintercept = median(tmrca$value), linetype = "dashed")


