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
metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025May22.csv')%>%
  filter(tipnames %in% nflg_mcc@phylo$tip.label) 

nflg_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/lineage_level_0_taxa/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_lineage0taxa_mcc.tree')

log_file <- beastio::readLog('./2025May22/global_analysis/global_subsampled_plain/lineage_level_0_taxa/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_lineage0taxa_1000.log',
                                 burnin = 0) %>%
  as_tibble() %>%
  select(starts_with('age.lineage')) %>%
  pivot_longer(everything(),
               values_to = 'draw',
               names_to = 'lineage') %>%
  mutate(lineage = gsub('age.lineage_level0_|\\.$', '', lineage))


################################### MAIN #######################################
# Map nodes between nextrain json (with autolin clades) and time-scaled BEAST tree
match <- as_tibble(ape::makeNodeLabel(lineage_json@phylo, method = "md5sum")) %>% 
  select(node, label) %>%
  rename(t1.node = node) %>%
  left_join(.,
            as_tibble(ape::makeNodeLabel(beast_mcc@phylo, method = "md5sum")) %>% 
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
  lapply(., function(x) as_tibble(beast_mcc)$node[as_tibble(beast_mcc)$label %in% x]) %>%
  lapply(., function(x) MRCA(beast_mcc, x)) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 1) %>%
  rename(lineage = name,
         root_node = value)

level_2_inferred <- level_2[names(level_2) %in% (all_levels_tbl %>% filter(is.na(root_node)) %>%  pull(lineage))] %>%
  lapply(., function(x) offspring(lineage_json, x , type = 'tips')) %>%
  lapply(., function(x) as_tibble(lineage_json)$label[as_tibble(lineage_json)$node %in% x]) %>%
  lapply(., function(x) as_tibble(beast_mcc)$node[as_tibble(beast_mcc)$label %in% x]) %>%
  lapply(., function(x) MRCA(beast_mcc, x)) %>%
  enframe() %>%
  unnest(value) %>%
  mutate(level = 2) %>%
  rename(lineage = name,
         root_node = value)

# Update all_lineage tibble with inferred root nodes (for missing only)
all_levels_tbl %<>% 
  rows_patch(level_1_inferred, by = c('lineage', 'level')) %>%
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



basic_tree <- nflg_mcc %>%
  left_join(metadata_in_tree %>% 
              dplyr::select(is_europe, tipnames) %>%
              rename(label = tipnames),
            by = 'label') %>%
  left_join(level_0_tbl, by = 'node') %>%
  ggtree(mrsd = '2024-10-12') + 
  theme_tree2() +
  scale_y_reverse(expand = expansion(mult = c(0, .05))) + 
  scale_x_continuous(limits = c(1900, 2031))

basic_tree + 
  #Backgrounds
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
  
  # Lineage TMRCA densities
  stat_slab(data = log_file , 
            aes(x =  draw,
                group = lineage,
                height = 250,
                slab_fill = lineage,
                slab_colour = lineage,
                y =550),
            slab_alpha = 0.7,
            #slab_colour = 'black',
            slab_linewidth = 0.75,
            #normalize = 'groups',
            p_limits = c(0.001, 0.999),
            inherit.aes = FALSE) +
  geom_linerange(data = log_file %>%
                 summarise(mu = mean(draw), .by = lineage) %>%
                 left_join(level_0_tbl) %>%
                 left_join(basic_tree[["data"]] %>% select(node, y)), 
               aes(x= mu,
                   ymin = Inf,
                   ymax = y,
                   group = lineage,
                   colour = lineage),
               inherit.aes = FALSE, linetype = "dashed", position = position_dodge(width = 0.5)) +
  geom_point2(aes(subset= !is.na(lineage), colour = lineage))  +
  scale_colour_manual(values = lineage_colours)+
  scale_colour_manual(aesthetics = 'slab_colour',values = lineage_colours)+
  
  scale_fill_manual(aesthetics = 'slab_fill',
                    values = lineage_colours)


write_csv(all_levels_tbl, './2025May22/nomenclature/mcc_lineage_root_nodes.csv')

#################################### END #######################################
################################################################################  