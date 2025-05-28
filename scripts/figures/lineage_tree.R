################################################################################
## Script Name:        Usutu Figure 1
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:      2025-05-27
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)

memory.limit(30000000)

############################## DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(tidytree)
library(TreeTools)
library(ggtree)
library(ggtreeExtra)
library(beastio)
library(ggmcmc)
library(ggnewscale)
################################### DATA #######################################
# Read and inspect data
nflg_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_mcc.tree')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025May22.csv')%>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

lineage_roots <- read_csv('./2025May22/nomenclature/mcc_lineage_root_nodes.csv')

all_node_heights <- read_delim('./2025May22/global_analysis/treestat_test') %>%
  pivot_longer(cols = starts_with('Height'), 
               values_to = 'draw',
               names_to = 'node') %>%
  mutate(node = str_extract(node, '\\d+') %>% as.double())

all_logs <- beastio::readTreeLog('./2025May22/global_analysis/global_subsampled_plain/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_1000.trees',
                                 burnin = 0)
lineage_tmrcas <- all_logs %>%
  lapply(., function(x) distRoot(x, 'PV523802|9731|homo_sapiens|HU|2024-10-12') -  distRoot(x, lineage_roots$root_node) %>%
           setNames(lineage_roots$root_node) %>% 
           as_tibble_row()) %>%
  
  bind_rows(.id = 'draw') %>%
  pivot_longer(-1, values_to = 'node_height',
               names_to = 'root_node') %>%
  mutate(root_node = as.numeric(root_node)) %>%
  left_join(lineage_roots)

################################### MAIN #######################################
# Main analysis or transformation steps
shading_intervals <- seq(1920, 2030, by = 10)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2030)

#a4243b, #be776b #d8c99b, #d8b06c, #d8973c, #cb7d36, #bd632f, #273e47

nflg_mcc %>%
  left_join(metadata_in_tree %>% 
              dplyr::select(is_europe, tipnames) %>%
              rename(label = tipnames),
            by = 'label') %>%
  ggtree(mrsd = '2024-10-12') + 
  theme_tree2() + 
  #scale_y_reverse() + 
  
  #Background
  annotate("rect", 
           xmin = shading_intervals[seq(1, length(shading_intervals), 2)], 
           xmax = shading_intervals_end[seq(1, length(shading_intervals_end), 2)], 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
  #Tip Colour (In/Out of Europe)
  geom_tippoint(aes(fill = as.factor(is_europe)),
                size = 3, 
                shape = 21) +
                  
  scale_fill_manual(values = c("0" = '#24a489', 
                                 '1' = '#2f91bd'), 
                      guide = 'none') +
  
  # Plot tmrca density
  new_scale_fill()+

  geom_density(data = lineage_tmrcas %>% filter(level == 0) %>% filter(lineage == 'A'), 
               aes(x =  2024.779-node_height,
                   fill = lineage,
                   y = after_stat(density)*1000),
              # bounds = GetNodeAges("2023-09-23",  c(70.5184216919521, 179.175463927162)) %>% decimal_date() %>% sort(),
               alpha = 0.5, 
               inherit.aes = FALSE) +
  
  geom_vline(xintercept = median(tmrca$value), linetype = "dashed")
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################  