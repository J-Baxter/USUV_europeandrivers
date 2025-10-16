################################################################################
## Script Name:        Extract posterior estimates from NFLG trees
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-08-15
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(treeio)
library(beastio)
library(ggmcmc)
library(ggdist)

################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("^NFLG", basename(dirs))]

logfilepaths <- sapply(dirs, 
                       list.files,
                       pattern = "(SG|constant)\\.log$",
                       full.names = TRUE, 
                       simplify = F) %>%
  Filter(length,.)%>%
  sapply(.,
         # If constant run is present, use this rather than skygrid
         function(x) ifelse(any(grepl('constant', x)), 
                            x[grepl('constant', x)], 
                            x),
         simplify = F)

tree_logs <- lapply(logfilepaths,
                    readLog,
                    burnin = 0.1) %>%
  lapply(., ggs) %>%
  bind_rows(., .id = 'clade') 


################################### MAIN #######################################
# Main analysis or transformation steps

# Summarise all reported parameters
median_hcdis <- tree_logs %>%
  group_by(clade, Parameter) %>%
  point_interval(.width = 0.95,
                 .point = median, 
                 .interval = hdci)

tmrca_summary <- median_hcdis %>%
  filter(grepl('^age.root', Parameter)) %>%
  dplyr::select(clade, starts_with('value')) %>%
  mutate(across(starts_with('value'), .fns = ~ date_decimal(.x))) %>%
  rename_with(~gsub('^value', 'TMRCA', .x) %>%
                gsub('\\.', '_', .),
              .cols = starts_with('value'))

median_hcdis %>%
  filter(grepl('^treeLength', Parameter)) %>%
  dplyr::select(clade, starts_with('value')) %>%
  #mutate(across(starts_with('value'), .fns = ~ date_decimal(.x))) %>%
  rename_with(~gsub('^value', 'TMRCA', .x) %>%
                gsub('\\.', '_', .),
              .cols = starts_with('value'))


rate_summary <- median_hcdis %>%
  filter(grepl('^meanRate$', Parameter)) %>%
  dplyr::select(clade, starts_with('value')) %>%
  rename_with(~gsub('^value', 'substitutionrate', .x) %>%
                gsub('\\.', '_', .),
              .cols = starts_with('value'))

summary_tbl <- left_join(tmrca_summary,
                         rate_summary)


################################### OUTPUT #####################################
# Save output files, plots, or results
write_csv(summary_tbl,
          './2025Jun24/europe_clusters/nflg_posterior_summary.csv')

#################################### END #######################################
################################################################################