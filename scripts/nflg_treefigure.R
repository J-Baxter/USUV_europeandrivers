####################################################################################################
####################################################################################################
## Script name: NFLG MCC USUV tree visualisation
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
library(ggtree)
library(ggnewscale)

# User functions


############################################## DATA ################################################
nflg_mcc <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_mcc.tree')
nflg_ca <- read.beast('./2024Oct20/test_beast/USUV_2024Oct20_nflg_subsample1_SRD06_RelaxLn_constant_ca.tree')

phyclip_lineages <- read_delim('./data/phyCLIP/cluster_optimal_parameter_cs3_fdr0.1_gam3.0_sol0_f0_zero-branch-length-collapsed_rooted_nflg.txt') %>%
  dplyr::select(c(1,2)) %>%
  rename(lineage = CLUSTER,
         label = TAXA) %>%
  mutate(lineage =as.factor(lineage))

metadata_in_tree <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv') %>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(nflg_mcc@phylo)) #sanity check


############################################## MAIN ################################################

most_recent_date <- metadata_in_tree %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym
  max(na.rm = TRUE)

ggtree(nflg_ca,
       mrsd = most_recent_date) + 
  theme_tree2() +

p <- nflg_mcc %>% 
  full_join(phyclip_lineages) %>%
  full_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames, lineage) %>%
              rename(label = tipnames),
            by = 'label') %>%
  ggtree(mrsd = most_recent_date) + 
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(shape = is.na(sequence_accession), fill = is.na(sequence_accession)),
                colour = NA) +
  scale_shape_manual(values = c("TRUE" = 23),'New Sequences') +
  scale_fill_manual(values = c("TRUE" = 'lightblue'), 'New Sequences') + 
  
  # node colour to show pp support
  geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  scale_color_distiller(palette = 'OrRd', direction = 1, 'Posterior Support') + 
  theme_tree2() +
  
  # add new scale for heat map
  new_scale_fill() 
  

# Add 'is europe'
p1 <- gheatmap(p + new_scale_fill(), metadata_in_tree %>% 
           dplyr::select(is_europe, tipnames) %>%
           mutate(is_europe = as.factor(is_europe)) %>%
           column_to_rownames(var = 'tipnames'), 
         offset=1, 
         width=0.05, 
         colnames = FALSE) +
  scale_fill_brewer('Location', palette = 'Dark2', labels = c('0' = 'Outside of Europe',
                                                        '1' = 'Within Europe')) +
  scale_x_ggtree(breaks = c(1920, 1940, 1960, 1980, 2000, 2020),
    labels = c(1920, 1940, 1960, 1980, 2000, 2020)) + 
  scale_y_continuous(expand=c(0, 0.3))

# Add phyclip lineages
p2 <- gheatmap(p1 + new_scale_fill(), phyclip_lineages %>% 
                 column_to_rownames(var = 'label'), 
               offset=8, 
               width=0.05, 
               colnames = FALSE,
               custom_column_labels = '') +
  #scale_fill_brewer('Lineage', palette = 'Paired') + #colour scheme needs factors in the right order so sublineages are the 'correct' colour
  scale_x_ggtree(breaks = c(1920, 1940, 1960, 1980, 2000, 2020),
                 labels = c(1920, 1940, 1960, 1980, 2000, 2020)) + 
  scale_y_continuous(expand=c(0, 0.3))

p3 <- gheatmap(p2 + new_scale_fill(), metadata %>% column_to_rownames(var = 'tipnames') %>% dplyr::select(lineage), 
               offset=16, 
               width=0.05, 
               colnames = FALSE,
               custom_column_labels = '') +
  #scale_fill_brewer('Lineage', palette = 'Paired') + #colour scheme needs factors in the right order so sublineages are the 'correct' colour
  scale_x_ggtree(breaks = c(1920, 1940, 1960, 1980, 2000, 2020),
                 labels = c(1920, 1940, 1960, 1980, 2000, 2020)) + 
  scale_y_continuous(expand=c(0, 0.3))
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################