####################################################################################################
####################################################################################################
## Script name: Genetic distance based clustering
##
## Purpose of script:
##
## Date created: 2024-10-28
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)
library(ade4)

# User functions


############################################## DATA ################################################
temp_tree <- read.tree('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta.treefile') %>%
  multi2di()

temp_tree$node.label <- ifelse(temp_tree$node.label == '', '100',temp_tree$node.label)

write.tree(temp_tree, './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1_bifur.tree')

aln_subsampled <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta',
                           format = 'fasta',
                           as.matrix = T)


metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')


metadata_in_tree <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv') %>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(nflg_mcc@phylo)) #sanity check

############################################## MAIN ################################################

europe_tips <- metadata_in_tree %>%
  filter(is_europe == 1) %>% 
  pull(tipnames) 


non_europe_tips <- metadata_in_tree %>%
  filter(is_europe == 0) %>% 
  pull(tipnames) 


############################################## WRITE ###############################################

geneticdistance_withineurope <- dist.dna(aln_subsampled[rownames(aln_subsampled) %in% europe_tips, ],
                                         model = "TN93")

geneticdistance_noneurope <- dist.dna(aln_subsampled[rownames(aln_subsampled) %in% non_europe_tips, ],
                                         model = "TN93")



within_europe <- geneticdistance_withineurope %>%
  as.matrix() %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  filter(!duplicated(paste0(pmax(a, b), pmin(a, b))))%>%
  mutate(type = 'Within Europe')


non_europe <- geneticdistance_noneurope %>%
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

all_genetic_distances <- bind_rows(non_europe,
                                     within_europe)


ggplot(all_genetic_distances ) +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
  scale_fill_brewer(NULL,
                    palette = 'Dark2')+
  scale_colour_brewer(NULL,
                      palette = 'Dark2')+
  scale_x_continuous('Genetic Distance',
                     expand = c(0.005,0))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic(base_size = 20)+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))
############################################## END #################################################
####################################################################################################
####################################################################################################
