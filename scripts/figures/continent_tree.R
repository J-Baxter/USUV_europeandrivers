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
library(treeio)
library(ggnewscale)
library(ggdist)
library(ggsci)
GetLineageRoots <- function(treedata, lineage_name){
  tree <- treedata@phylo
  data <- as_tibble(treedata) %>% 
    select(label, {{lineage_name}})
  
  # Infer all subtrees and attach associated metadata
  out <- subtrees(tree) %>% 
    lapply(., function(x) as.treedata(x) %>%
             left_join(data,
                       by= 'label')) %>%
    
    # Convert to tibble format
    lapply(., as_tibble) %>%
    bind_rows(.id = 'subtree') %>%
    group_by(subtree) %>%
    
    # Include only subtrees that are monophyletic wrt lineage
    filter(n_distinct(.data[[lineage_name]]) == 1)  %>% 
    
    # Exclude unassigned subtrees
    filter(.data[[lineage_name]] != 'not assigned') %>%
    
    mutate(n = n()) %>%
    ungroup() %>%
    group_by(.data[[lineage_name]]) %>%
    # include only the largest subtree
    filter(n == max(n)) %>%
    
    # extract root node as list
    filter(node == parent) %>%
    
    select(label,
           {{lineage_name}}) %>%
    summarise(named_vec = list(label), .groups = "drop") %>%
    deframe()
  
  return(out)
}



################################### DATA #######################################
# Read and inspect data
#nflg_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/lineage_level_0_taxa/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_lineage0taxa_mcc.tree')
nflg_hipstr <- read.beast('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_constant_hipstr.tree')


#og_file <- beastio::readLog('./2025May22/global_analysis/global_subsampled_plain/lineage_level_0_taxa/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_lineage0taxa_1000.log',
                             #burnin = 0) %>%
  #as_tibble() %>%
  #select(starts_with('age.lineage')) %>%
  #pivot_longer(everything(),
               #values_to = 'draw',
               #names_to = 'lineage') %>%
  #mutate(lineage = gsub('age.lineage_level0_|\\.$', '', lineage))


log_file <- beastio::readLog('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_constant.log',
                             burnin = 0) %>%
  as_tibble() %>%
  select(starts_with('age.lineage')) %>%
  pivot_longer(everything(),
               values_to = 'draw',
               names_to = 'lineage') %>%
  mutate(lineage = gsub('age.lineage_level0_|\\.$', '', lineage))


metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')%>%
  filter(tipnames %in% nflg_hipstr@phylo$tip.label) 

lineage_json <- read.nextstrain.json('./2025Jun24//nomenclature/subsample/USUV_2025Jun24_lineages.json')

lineage_tbl <- read_csv('./2025Jun24/nomenclature/wi') %>%
  mutate(across(starts_with('GRI'), .fns = ~ gsub('not assigned', NA_character_, .x)))

root_lineage <- read_csv('./2025May22/nomenclature/mcc_lineage_root_nodes.csv')

################################### MAIN #######################################
# Map nodes between nextrain json (with autolin clades) and time-scaled BEAST tree
match <- as_tibble(ape::makeNodeLabel(lineage_json@phylo, method = "md5sum")) %>% 
  select(node, label) %>%
  rename(t1.node = node) %>%
  left_join(.,
            as_tibble(ape::makeNodeLabel(nflg_hipstr@phylo, method = "md5sum")) %>% 
              select(node,label) %>% 
              rename(t2.node=node)) %>%
  left_join(as_tibble(lineage_json) %>%
              select(t1.node = node, 
                     starts_with('GRI')))


# Infer Root nodes for lineage levels
level_0 <- GetLineageRoots(lineage_json, 'GRI Lineage Level 0')
level_1 <- GetLineageRoots(lineage_json, 'GRI Lineage Level 1') %>%
  { .[grepl("\\.", names(.))] }
level_2 <- GetLineageRoots(lineage_json, 'GRI Lineage Level 2') %>%
  { .[grepl("\\.\\d+\\.", names(.))] }


# Extract BEAST MCC root nodes for each lineage
level_0_tbl <- lapply(level_0, function(x) as_tibble(lineage_json)$node[as_tibble(lineage_json)$label == x]) %>%
  lapply(., function(x) match$t2.node[match$t1.node == x]) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 0) %>%
  rename(lineage = name,
         root_node = value)

level_1_tbl <- lapply(level_1, function(x) as_tibble(lineage_json)$node[as_tibble(lineage_json)$label == x]) %>%
  lapply(., function(x) match$t2.node[match$t1.node == x]) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 1) %>%
  rename(lineage = name,
         root_node = value)

level_2_tbl <- lapply(level_2, function(x) as_tibble(lineage_json)$node[as_tibble(lineage_json)$label == x]) %>%
  lapply(., function(x) match$t2.node[match$t1.node == x]) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 2) %>%
  rename(lineage = name,
         root_node = value)

all_levels_tbl <- bind_rows(level_0_tbl,
                            level_1_tbl,
                            level_2_tbl)


# Infer missing (no direct match between nodes)
level_1_inferred <- level_1[names(level_1) %in% (all_levels_tbl %>% filter(is.na(root_node)) %>%  pull(lineage))] %>%
  lapply(., function(x) offspring(lineage_json, x , type = 'tips')) %>%
  lapply(., function(x) as_tibble(lineage_json)$label[as_tibble(lineage_json)$node %in% x]) %>%
  lapply(., function(x) as_tibble(nflg_hipstr)$node[as_tibble(nflg_hipstr)$label %in% x]) %>%
  lapply(., function(x) MRCA(nflg_hipstr, x)) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 1) %>%
  rename(lineage = name,
         root_node = value)

level_2_inferred <- level_2[names(level_2) %in% (all_levels_tbl %>% filter(is.na(root_node)) %>%  pull(lineage))] %>%
  lapply(., function(x) offspring(lineage_json, x , type = 'tips')) %>%
  lapply(., function(x) as_tibble(lineage_json)$label[as_tibble(lineage_json)$node %in% x]) %>%
  lapply(., function(x) as_tibble(nflg_hipstr)$node[as_tibble(nflg_hipstr)$label %in% x]) %>%
  lapply(., function(x) MRCA(nflg_hipstr, x)) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 2) %>%
  rename(lineage = name,
         root_node = value)

# Update all_lineage tibble with inferred root nodes (for missing only)
all_levels_tbl %<>% 
  rows_patch(level_1_inferred, by = c('lineage', 'level')) %>%
  rows_patch(level_2_inferred, by = c('lineage', 'level'))


level_1_tbl %<>% 
  rows_patch(level_1_inferred, by = c('lineage', 'level')) %>%
  mutate(root_node = if_else(lineage == 'B.2', 977, root_node))

level_2_tbl  %<>% 
  rows_patch(level_2_inferred, by = c('lineage', 'level'))

################################### OUTPUT #####################################
# Plot on Tree
shading_intervals <- seq(1900, 2030, by = 10)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2030)
lineage_colours <- c('A' = '#a4243b',
                     'B' = '#be776b',
                     'C' = '#d8c99b',
                     'D' = '#d8b06c', 
                     # 'E' =  '#d8973c', 
                     'E' = '#cb7d36', 
                     'F' = '#bd632f', 
                     'G' = '#273e47')



nflg_hipstr %>%
  left_join(metadata_in_tree %>% 
              dplyr::select(is_europe, tipnames, nuts0_id, host_class) %>%
              rename(label = tipnames),
            by = 'label') %>% 
  left_join(level_0_tbl, by = join_by(node ==root_node)) %>% 
  left_join(lineage_tbl, by = 'label') %>% 
  ggtree(mrsd = '2024-10-12') + 
  theme_tree2() +
  #scale_y_reverse(expand = expansion(mult = c(0.02, .02))) +
  scale_x_continuous(breaks = seq(1905, 2025, by = 10)) + 
  #Backgrounds
  annotate("rect", 
           xmin = shading_intervals[seq(1, length(shading_intervals), 2)], 
           xmax = shading_intervals_end[seq(1, length(shading_intervals_end), 2)], 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
 # Tip Colour (In/Out of Europe)
  geom_tippoint(aes(fill = factor(is_europe)),
  size = 3, 
  shape = 21) +
  
  scale_fill_manual( values = c("0" = '#24a489', 
   '1' = '#2f91bd'), 
   labels = c('0' = 'Non-Europe',
              '1' = 'Europe'),
  lineage_colours) +
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = nuts0_id),
             width = 4,
             colour = "white",
             pwidth = 1.2,
             offset = 00.04) +
  scale_fill_d3(name = 'Country', 
                palette ='category20', 
                alpha = 0.99, 
                labels = c('AT' = 'Austria',
                           'BE' = 'Belgium',
                           'CF' = 'Central African Republic',
                           'DE' = 'Germany',
                           'ES' = 'Spain',
                           'FR' = 'France',
                           'HU' = 'Hungary',
                           'IL' = 'Israel',
                           'IT' = 'Italy',
                           'LU' = 'Luxembourg',
                           'NL' = 'Netherlands',
                           'PT' = 'Portugal',
                           'RS' = 'Serbia',
                           'SE' = 'Sweden',
                           'SK' = 'Slovakia',
                           'SN' = 'Senegal',
                           'UG' = 'Uganda',
                           'UK' = 'United Kingdom',
                           'EL' = 'Greece',
                           'ZA' = 'South Africa' ),
                guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 2, order = 1)) +
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = host_class),
             width = 4,
             colour = "white",
             pwidth = 1.2,
             offset = 00.04) +
  scale_fill_brewer('Host Class', 
                    labels = c('Aves' = 'Bird',
                               'insecta'= 'Insect',
                               'mammalia' = 'Mammal',
                               'arachnida' = 'Arachnid'),
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 3))+
  theme(legend.position = c(0.2,0.6),
        # legend.position = "bottom",       # Place legends at the bottom
        legend.box = "vertical",
        legend.direction = 'vertical'
  ) +

  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position.inside = c(0,1),
        legend.justification=c(0,1),
        text = element_text(size = 16)) 

                    
ggsave('~/Downloads/europe_tree.pdf', height = 30, width = 20, units = 'cm', dpi =360, device = "pdf")

write_csv(all_levels_tbl, './2025May22/nomenclature/mcc_lineage_root_nodes.csv')



#################################### END #######################################
################################################################################  