################################################################################
## Script Name:        USUV Nomenclature (May 2025)
## Purpose:            1) reconcile labelled json and BEAST MCC
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


################################### DATA #######################################
# Read and inspect data
lineage_json <- treeio::read.nextstrain.json('./2025Jun24/nomenclature/USUV_2025Jun24_lineages.json')


################################### MAIN #######################################
# Main analysis or transformation steps
taxa_for_beast <- as_tibble(lineage_json) %>%
  select(label, `GRI Lineage Level 0`) %>% 
  rename(taxa = `GRI Lineage Level 0`) %>%
  mutate(taxa = if_else(taxa == 'not assigned', NA_character_, taxa)) %>%
  drop_na(taxa) %>%
  filter(!grepl('^NODE', label)) %>%
  group_split(taxa) %>%
  set_names(LETTERS[1:length(.)]) %>%
  lapply(., function(x) x %>% select(label))


################################### OUTPUT #####################################
# Save output files, plots, or results
write_delim(as_tibble(lineage_json) %>%
            select(label, starts_with('GRI')) %>% 
            filter(!grepl('^Node', label)) ,
          './2025Jun24/nomenclature/wide_lineages.txt')


taxa_for_beast %>%
  mapply(function(x,y) write_delim(x, paste0('./2025Jun24/nomenclature/lineage_level0_', y, '.txt'),col_names = FALSE,),
         .,
         LETTERS[1:length(.)],
         SIMPLIFY = FALSE)

#################################### END #######################################
################################################################################  