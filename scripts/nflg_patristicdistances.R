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

# User functions


############################################## DATA ################################################
nflg_ca <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_ca.tree')
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')

phyclip_lineages <- read_delim('./data/phyCLIP/cluster_optimal_parameter_cs3_fdr0.1_gam3.0_sol0_f0_zero-branch-length-collapsed_rooted_nflg.txt') %>%
  dplyr::select(c(1,2)) %>%
  rename(lineage = CLUSTER,
         label = TAXA) %>%
  mutate(lineage =as.factor(lineage))


nflg_ca %<>% 
  left_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames) %>%
              rename(label = tipnames),
            by = 'label') 


phyclip_lineages <- read_delim('./data/phyCLIP/.txt') %>%
  dplyr::select(c(1,2)) %>%
  rename(lineage = CLUSTER,
         label = TAXA) %>%
  mutate(lineage =as.factor(lineage))

############################################## MAIN ################################################

europe_tips <- nflg_ca %>%
  as_tibble(nflg_ca) %>% 
  filter(is_europe == 1) %>% 
  pull(label) 
  

non_europe_tips <- nflg_ca %>%
  as_tibble(nflg_ca) %>% 
  filter(is_europe == 0) %>% 
  pull(label) 


patristic_distances_withineurope <- distTips(nflg_ca@phylo,
                                             tips = europe_tips,
                                             useC = FALSE) 


patristic_distances_noneurope <- adephylo::distTips(nflg_ca@phylo,
                                          tips = non_europe_tips,
                                          method = "patristic",
                                          useC = FALSE) 


within_europe <- patristic_distances_withineurope %>%
  as.matrix() %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b))))%>%
  mutate(type = 'Within Europe')


non_europe <- patristic_distances_noneurope %>%
  as.matrix()%>%
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
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  tipnames) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b)))) 

all_patristic_distances <- bind_rows(non_europe,
                                     within_europe)

ggplot(all_patristic_distances ) +
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
 

   
# Europe tips only

minimum_import_clades <- offspring(nflg_mcc, c(737, 624, 450, 852, 435), tiponly = TRUE) %>%
  lapply(., function(x) dplyr::select(as_tibble(nflg_mcc), c(label, node)) %>% 
           filter(node %in% x) %>%
           pull(label)) %>%
  set_names(LETTERS[1:5]) %>%
  lapply(., as_tibble) %>%
  bind_rows(,.id = 'cluster') %>%
  rename(tipnames = value) %>%
  left_join(metadata %>% dplyr::select(tipnames, is_europe)) %>%
  filter(is_europe == 1 ) %>%
  dplyr::select(-is_europe) %>%
  rename(label = tipnames)


within_europe_clades <- patristic_distances_withineurope %>%
  as.matrix() %>%
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



non_europe <- patristic_distances_noneurope %>%
  as.matrix() %>%
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
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  tipnames) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b)))) 

all_patristic_distances_clade <- bind_rows(non_europe,
                                     within_europe_clades)

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

within_europe_clades %>%
  ggplot() +
  geom_tile(aes(x = label_a, y = label_b, fill = distance)) + 
  scale_fill_viridis_c(option = 'C')

all_patristic_distances_clade %>%
  summarise(dist = m(distance), .by = c(lineage_a, lineage_b, type)) %>%
  drop_na(starts_with('lineage')) %>%
  filter(grepl('Between Clade', type)) %>%
  ggplot() +
  geom_tile(aes(x = lineage_a, y = lineage_b, fill = dist)) + 
  scale_fill_viridis_c(option = 'C')


hclust_median <- hclust(patristic_distances_withineurope, method = 'complete')
cut_avg <- list(cbind.data.frame(cluster= cutree(hclust_median, k = 5),  k = 5),
                     cbind.data.frame(cluster= cutree(hclust_median, k = 6),  k = 6),
                     cbind.data.frame(cluster= cutree(hclust_median, k = 7),  k = 7),
                     cbind.data.frame(cluster= cutree(hclust_median, k = 8),  k = 8),
                     cbind.data.frame(cluster= cutree(hclust_median, k = 9), k = 9),
                     cbind.data.frame(cluster= cutree(hclust_median, k = 10),  k = 10)) %>%
  lapply(., as_tibble, rownames = 'label') %>%
  bind_rows()
  

within_europe_clades <- patristic_distances_withineurope %>%
  as.matrix() %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(cut_avg %>% filter(k==8)) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(cluster) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, cluster)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('cluster'))

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


nflg_mcc %>% 
  full_join(cut_avg %>% filter(k==6),
            by = 'label') %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = as.factor(cluster)))
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################