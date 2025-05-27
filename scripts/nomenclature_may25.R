################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
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
library(treeio)
library(TreeTools)
library(ggtree)
library(ggtreeExtra)

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
lineage_json <- treeio::read.nextstrain.json('./2025May22/nomenclature/usuv_lineages.json')

autolin_labels <- read_tsv('./2025May22/global_analysis/nextstrain/automated-lineage-json/labels.tsv') %>%
  rename(label = sample)

beast_mcc <- read.beast('./2025May22/global_analysis/global_subsampled_plain/USUV_2025May22_noFLI_NFLG_subsampled_SRD06_RelaxLn_constant_mcc.tree')


################################### MAIN #######################################
# Main analysis or transformation steps
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
# Save output files, plots, or results

#################################### END #######################################
################################################################################  