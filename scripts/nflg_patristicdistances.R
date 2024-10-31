####################################################################################################
####################################################################################################
## Script name: NFLG patristic distances between/within Europe
##
## Purpose of script:
##
## Date created: 2024-10-25
##
##
########################################## SYSTEM OPTIONS ##########################################


  
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
library(ggtreeExtra)

# User functions


############################################## DATA ################################################
nflg_ca <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_ca.tree')
nflg_mcc <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_mcc.tree')
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')
patristic_distance_clusters <- read_csv('./2024Oct20/nflg_patristic_distance_clusters.csv')


nflg_ca %<>% 
  left_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames) %>%
              rename(label = tipnames),
            by = 'label') 


############################################## MAIN ################################################

# Calculate pairwise patristic distances for all tips on phylogeny
patristic_distances <- distTips(nflg_ca@phylo,
                                method = "patristic") 

# vector of Europe tips
europe_tips <- nflg_ca %>%
  as_tibble(nflg_ca) %>% 
  filter(is_europe == 1) %>% 
  pull(label) 
  

# vector of non-Europe tips
non_europe_tips <- nflg_ca %>%
  as_tibble(nflg_ca) %>% 
  filter(is_europe == 0) %>% 
  pull(label) 


# 1: Compare within-europe and outside eureop
within_europe <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b))))%>%
  mutate(type = 'Within Europe')


outside_europe_within <- patristic_distances %>%
  as.matrix() %>%
  .[non_europe_tips, non_europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'tipnames') %>%
  left_join(metadata %>% dplyr::select(tipnames, nuts0_id)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(nuts0_id) == 1 ~ 'Outside Europe: Within Country',
                              .default = 'Outside Europe: Between Country')) %>%
  ungroup() %>%
  filter(type == 'Outside Europe: Within Country') %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  tipnames) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b)))) 


outside_europe_between <- patristic_distances %>%
  as.matrix() %>%
  .[non_europe_tips, non_europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'tipnames') %>%
  left_join(metadata %>% dplyr::select(tipnames, nuts0_id)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(nuts0_id) == 1 ~ 'Outside Europe: Within Country',
                          .default = 'Outside Europe: Between Country')) %>%
  ungroup() %>%
  filter(type == 'Outside Europe: Between Country') %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  tipnames) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b)))) 




 p1b <-  bind_rows(within_europe,
            outside_europe_within,
            outside_europe_between) %>%
  ggplot() +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
  scale_fill_brewer(NULL,
                    palette = 'Dark2')+
  scale_colour_brewer(NULL,
                      palette = 'Dark2')+
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic(base_size = 20)+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))
 
nflg_mcc <- read.beast('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits_mcc.tree')
most_recent_date <- metadata_in_tree %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym
  max(na.rm = TRUE)


# Define the shading intervals - age root for tmrca. Not sure yet on other nodes
shading_intervals <- seq(1920, 2030, by = 10)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2030)

p1 <- nflg_mcc %>% 
  ggtree(mrsd = most_recent_date, aes(colour = is_europe)) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20)) +
  scale_colour_brewer('Location', palette = 'Dark2', labels = c('not_europe' = 'Outside of Europe',
                                                                'europe' = 'Within Europe'))+
  new_scale_colour()+
  geom_nodepoint(aes(colour = is_europe))+
  scale_colour_brewer('Location', palette = 'Dark2', labels = c('not_europe' = 'Outside of Europe',
                                                                'europe' = 'Within Europe')) + 
  #Background
  annotate("rect", 
           xmin = shading_intervals[seq(1, length(shading_intervals), 2)], 
           xmax = shading_intervals_end[seq(1, length(shading_intervals_end), 2)], 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") 

cowplot::plot_grid( p1, p1b, align = 'hv')

# Distributions for minimum number of imports (ie, one intro per phyly of EU sequences)
minimum_import_clades <- offspring(nflg_mcc, c(737, # Africa 3
                                               449, #EU1 + 2 
                                               852, #Africa 2
                                               435 #Africa 2 - KU760915
                                               ), tiponly = TRUE) %>%
  lapply(., function(x) dplyr::select(as_tibble(nflg_mcc), c(label, node)) %>% 
           filter(node %in% x) %>%
           pull(label)) %>%
  set_names(LETTERS[1:4]) %>%
  lapply(., as_tibble) %>%
  bind_rows(.,.id = 'cluster') %>%
  rename(tipnames = value) %>%
  left_join(metadata %>% dplyr::select(tipnames, is_europe)) %>%
  filter(is_europe == 1 ) %>%
  dplyr::select(-is_europe) %>%
  rename(label = tipnames)


within_europe_mic <-  patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(minimum_import_clades) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(cluster) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, cluster)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('cluster'))



bind_rows(within_europe_mic,
          outside_europe_within,
          outside_europe_between) %>%
  ggplot() +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
  scale_fill_brewer(NULL,
                    palette = 'Dark2')+
  scale_colour_brewer(NULL,
                      palette = 'Dark2')+
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic(base_size = 20)+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))




# Using patristic distance clades from patristic_distance_clusters.R
pdist30_clades <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(patristic_distance_clusters %>% dplyr::select(label, dist_30)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(dist_30) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, dist_30)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('dist_30'))


pdist40_clades <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(patristic_distance_clusters %>% dplyr::select(label, dist_40)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(dist_40) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, dist_40)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('dist_40'))


pdist50_clades <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(patristic_distance_clusters %>% dplyr::select(label, dist_50)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(dist_50) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, dist_50)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('dist_50'))



pdist60_clades <- patristic_distances %>%
    as.matrix() %>%
    .[europe_tips,europe_tips] %>%
    as_tibble(rownames = 'a') %>%
    pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
    filter(a != b) %>%
    rowid_to_column(var = 'pair_id') %>%
    pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
    left_join(patristic_distance_clusters %>% dplyr::select(label, dist_60)) %>%
    group_by(pair_id) %>%
    mutate(type = case_when(n_distinct(dist_60) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, dist_60)) %>%
    filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
    drop_na(starts_with('dist_60'))
    
  
  

ggplot(within_europe_clades) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
  scale_fill_brewer(NULL,
                    palette = 'Dark2')+
  scale_colour_brewer(NULL,
                      palette = 'Dark2')+
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic(base_size = 20)+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))


############################################## WRITE ###############################################
# Panel 1: left column annotated phylo, right column corresponding patristic distance distribution


# Phylos
phylo_60 <- nflg_mcc %>% 
  left_join(patristic_distance_clusters) %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_60)) +
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL,
                  na.value = NA) + 
  theme_tree2() +
  theme(legend.position = 'none')

phylo_50 <- nflg_mcc %>% 
  left_join(patristic_distance_clusters) %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_50)) +
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL,
                  na.value = NA) + 
  theme_tree2() +
  theme(legend.position = 'none')

phylo_40 <- nflg_mcc %>% 
  left_join(patristic_distance_clusters) %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_40)) +
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL,
                  na.value = NA) + 
  theme_tree2() +
  theme(legend.position = 'none')

mycolours <- c("#1B9E77" ,"#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#4f77a4" ,"#48aa35","#666666")
phylo_30 <- nflg_mcc %>% 
  left_join(patristic_distance_clusters) %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_30)) +
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL,
                  na.value = NA) + 
  theme_tree2() +
  theme(legend.position = 'none')

phylo_grid <- plot_grid(phylo_60,
                        phylo_50,
                        phylo_40,
                        phylo_30,
                        align = 'hv',
                        ncol = 4)


# Distributions

dist_60 <- ggplot(pdist60_clades) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
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

dist_40 <- ggplot(pdist40_clades) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
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
  theme(legend.position = 'none')

dist_50 <- ggplot(pdist50_clades) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
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
  theme(legend.position = 'none')

dist_30 <- ggplot(pdist30_clades) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
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
  theme(legend.position = 'none')


dist_grid <- plot_grid(dist_60,
                        dist_50,
                        dist_40,
                        dist_30,
                        align = 'hv',
                        ncol= 4)

plot_grid(phylo_grid,
          dist_grid,
          #ncol = 4,
          nrow = 2,
          align = 'hv')
############################################## END #################################################
####################################################################################################
####################################################################################################