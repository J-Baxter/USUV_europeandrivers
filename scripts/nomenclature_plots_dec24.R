####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-05-26
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(treeio)
library(TreeTools)
library(ggtree)
library(ggtreeExtra)
# User functions


############################################## DATA ################################################
test <- treeio::read.nextstrain.json('./2025May22/nomenclature/usuv_lineages.json')
autolin_labels <- read_tsv('./2025May22/global_analysis/nextstrain/automated-lineage-json/labels.tsv') %>%
  rename(label = sample)
beast_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_mcc.tree')

plain_ml <- ape::read.tree('./2025May22/global_analysis/global_subsampled_plain/USUV_2025May22_alldata_aligned_formatted_noFLI_NFLG_ml.tree') %>%
  midpoint_root()


as_tibble(test) %>%
  mutate(ancestors = )

t1 <- test %>% 
  as_tibble() %>% 
  select(c(1,2,4))

t2 <- beast_mcc %>% 
  as_tibble() %>% 
  select(c(1,2,4)) %>%
  rename(new_parent = parent,
         new_node = node)

# match tips
terminal <- t1 %>% filter(node<494) %>%
  left_join(t2 %>% drop_na,  
              by = 'label') 
  

# match parents of tips
t4 <- t1 %>% filter(node>=494) %>%
  left_join(terminal %>% select(node = parent, new_node = new_parent)) %>%
  left_join(t2 %>% select(-label), by = join_by(new_node)) %>%
  drop_na()

# match parents of parents

t1 %>% filter(node>=494) %>%
  left_join(terminal %>% select(node = parent, new_node = new_parent)) %>%
  left_join(t2 %>% select(-label), by = join_by(new_node)) %>%
  filter(is.na(new_node)) %>%
  select(-starts_with('new')) %>%
  left_join(t4 %>% select(node = parent,
                          new_node = new_parent), by = join_by(node)) %>%
  left_join(t2 %>% select(-label), by = join_by(new_node)) %>%
  
  
  left_join(t2 %>% select(node, new_parent = parent),  
            by = join_by(new_node==node)) %>%
  drop_na()
  
PathtoTips <- function(tree, node, tip){
  ancestors <- ancestor(tree, tip)
  branch_start <- tree$edge[1,]
  branch_length <- tree$edge.length
  
  dist_to_tip <- 0
  
  for (anc in ancestors){
    while(anc != node){
      index = which(branch_start==anc)
      dist_to_tip = dist_to_tip + branch_length[index]
    }
  }
  return(dist_to_tip)
}
ancestors <- ancestor(plain_ml, 1)
branch_start <- plain_ml$edge[,1]
branch_length <- plain_ml$edge.length

dist_to_tip <- 0

  for (anc in ancestors){
  print(anc)
  index = which(branch_start==anc)
  print(index)
  dist_to_tip = dist_to_tip + branch_length[index]
  
  }

PathtoTips(plain_ml, 702, 1)
as_tibble(plain_ml) %>%
  mutate(N = length(offspring(plain_ml, node, tiponly = TRUE)),
         D = branch.length,
         S = )


############################################## MAIN ################################################
primary_lineages <- autolin_labels %>%
  dplyr::select(lineage) %>%
  filter(!grepl('\\.', lineage)) %>%
  distinct() %>%
  rename(level_0 = lineage)

test_mrca_level0 <- test@data %>%
  filter(`GRI Lineage Level 0` != 'not assigned') %>%
  group_by(`GRI Lineage Level 0`) %>%
  slice_min(num_date) %>% 
  dplyr::select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 0`)


test_mrca_level1 <- test@data %>%
  filter(`GRI Lineage Level 1` != 'not assigned') %>%
  group_by(`GRI Lineage Level 1`) %>%
  slice_min(num_date) %>% 
  dplyr::select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 1`) %>%
  mutate(label = case_when(!grepl('\\.', label) ~ NA_character_, .default = label)) %>%
  drop_na()


test_mrca_level2 <- test@data %>%
  filter(`GRI Lineage Level 2` != 'not assigned') %>%
  group_by(`GRI Lineage Level 2`) %>%
  slice_min(num_date) %>% 
  dplyr::select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 2`) %>%
  mutate(label = case_when(!grepl('\\.[[:digit:]]\\.', label) ~ NA_character_, .default = label)) %>%
  drop_na()

old_nomenclature <- as_tibble(cbind.data.frame(
  label = c('EU1', 'EU2', 'EU3', 'EU4', 'EU5', 'AF1', 'AF2', 'AF3'),
  node = c(688, 691, 621, 619, 604, 485, 480, 496))) %>%
  drop_na()



# set colour scheme
level_0_cols <- pal_d3(palette = "category10")(nrow(primary_lineages))

primary_lineages %<>%
  cbind.data.frame(level_0_cols)

level_1_cols <- autolin_labels %>%
  dplyr::select(lineage) %>%
  distinct() %>%
  #filter(grepl('\\.', lineage)) %>%
  separate_wider_delim(lineage, delim = '.', names = c('level_0', 'level_1', 'level_2'),
                       too_few = 'align_start') %>%
  dplyr::select(-level_2) %>%
  distinct() %>%
  mutate(level_1  = as.integer(level_1)) %>%
  group_by(level_0) %>%
  arrange(level_1,.by_group = TRUE) %>%
  filter(!(is.na(level_1) & n() > 2)) %>%
  mutate(n = cur_group_id()) %>%
  mutate(level_1_col = lighten(level_0_cols[n], 
                               amount = seq(0.2, 0.8, length.out = n()))) %>%
  ungroup() %>%
  drop_na() %>%
  pull(level_1_col)


level_2_cols <- autolin_labels %>%
  dplyr::select(lineage) %>%
  distinct() %>%
  #filter(grepl('\\.', lineage)) %>%
  separate_wider_delim(lineage, delim = '.', names = c('level_0', 'level_1', 'level_2'),
                       too_few = 'align_start') %>%
  dplyr::select(-level_1) %>%
  distinct() %>%
  mutate(level_2  = as.integer(level_2)) %>%
  group_by(level_0) %>%
  arrange(level_2,.by_group = TRUE) %>%
  filter(!(is.na(level_2) & n() > 2)) %>%
  mutate(n = cur_group_id()) %>%
  mutate(level_2_col = lighten(level_0_cols[n], 
                               amount = seq(0.2, 0.8, length.out = n()))) %>%
  ungroup() %>%
  drop_na() %>%
  pull(level_2_col)



# plot NFLG nextstrain (not time scaled) phylogeny
p1 <- test %>% 
  ggtree() + 
  geom_cladelab(data = test_mrca_level0 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                textcolour = NA,
                barsize = 15,
                offset = 0.004) + 
  
  geom_cladelab(data = test_mrca_level0 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 0,
                offset = 0.004,
                offset.text = -0.0004) + 
  scale_colour_manual(values = level_0_cols) + 
  new_scale_colour()+
  
  geom_cladelab(data = test_mrca_level1 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                textcolour = NA,
                barsize = 15,
                offset = 0.0025) + 
  
  geom_cladelab(data = test_mrca_level1 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 0,
                offset = 0.0025,
                offset.text = -0.00055) + 
  scale_colour_manual(values = level_1_cols) + 
  
  new_scale_colour() +
  geom_cladelab(data = test_mrca_level2 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 0) + 
  
  theme(legend.position = 'none') +
  new_scale_colour() +
  geom_tiplab(align = TRUE, aes(label = "", colour =`GRI Lineage Level 0` )) +
  scale_colour_manual(values = c(level_0_cols, 'white')) + 
  
  # old
  new_scale_colour() +
  geom_cladelab(data = old_nomenclature, 
                mapping = aes(node = node, label = label, colour = label, alpha = 0.05), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = NA,
                barsize = 15,
                offset = 0.006) + 
  scale_colour_brewer(palette = 'Pastel2')+
  
  geom_cladelab(data = old_nomenclature , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                offset = 0.006,
                textcolour = 'black',
                offset.text = -0.0006,
                barsize = 0) 


# collapse level 2 lineages
p <-purrr::reduce(
  test_mrca_level2$node,
  \(x,y) 
  collapse(x, y, mode = "max", fill = "black", size = 0.05),
  .init = p1
)


############################################## WRITE ###############################################
ggsave(p, filename = "~/Downloads/nomenclature_05.png", width = 13, height = 25)


############################################## END #################################################
####################################################################################################
####################################################################################################