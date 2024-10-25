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
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(treeio)
library(TreeTools)
library(ape)
library(ad4)
library(adephylo)

# User functions


############################################## DATA ################################################
nflg_ca <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_ca.tree')

nflg_ca %<>% 
  full_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames) %>%
              rename(label = tipnames)) %>%
  
  mutate(branch)
  
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
                                             useC = FALSE) %>%
  as.matrix()

patristic_distances_noneurope <- adephylo::distTips(nflg_ca@phylo,
                                          tips = non_europe_tips,
                                          method = "patristic",
                                          useC = FALSE) %>%
  as.matrix()


within_europe <- patristic_distances_withineurope %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b))))%>%
  mutate(type = 'Within Europe')


non_europe <- patristic_distances_noneurope %>%
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
within_europe_clades <- patristic_distances_withineurope %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
  left_join(phyclip_lineages) %>%
  group_by(pair_id) %>%
  mutate(type = case_when(n_distinct(lineage) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, type), names_from = 'name', values_from =  c(label, lineage)) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b)))) %>%
  drop_na(starts_with('lineage'))



non_europe <- patristic_distances_noneurope %>%
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

ggplot(all_patristic_distances_clade) +
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

all_patristic_distances_clade %>%
  summarise(dist = m(distance), .by = c(lineage_a, lineage_b, type)) %>%
  drop_na(starts_with('lineage')) %>%
  filter(grepl('Between Clade', type)) %>%
  ggplot() +
  geom_tile(aes(x = lineage_a, y = lineage_b, fill = dist)) + 
  scale_fill_viridis_c(option = 'C')

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################