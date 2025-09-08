################################################################################
## Script Name:       Prior distributions for NFLG+partial scripts
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-08-19
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
library(ape)

formdate <- function(x){
  if(nchar(x) == 4){
    out <- ymd(x, truncated = 2)
  }else if(nchar(x) == 7){
    out <- ymd(x, truncated = 1)
  }else{
    out <- ymd(x)
  }
  
  return(decimal_date(out))
}


################################### DATA #######################################
# Read and inspect data
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("^NFLG", basename(dirs))]

logfilepaths <- sapply(dirs, 
                       list.files,
                       pattern = "(SG|constant|_test)\\.log$",
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


temp <- list.files('./2025Jun24/alignments',
                   pattern = 'partial_[:A-Z:]{1,4}_subsampled_filtered',
                   full.names = T) 


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
                         rate_summary)  %>%
  mutate(across(starts_with('TMRCA'), .fns = ~decimal_date(.x))) %>% 
  mutate(clade = str_split_i(clade, '_', 3)) 


temp %>%
  lapply(read.dna, format = 'fasta', as.matrix = T) %>%
  lapply(., rownames) %>%
  lapply(., function(x) str_extract(x, '\\d{4}(-\\d{2}){0,1}(-\\d{2}){0,1}$')) %>%
  lapply(., function(x) sapply(x, formdate, simplify = F) %>% unname(.) %>% unlist(.)) %>%
  lapply(., function(x) x[which.max(x)]) %>% 
  setNames(as.roman(2:(length(.)+1))) %>%
  enframe(name = 'clade',
          value = 'mrd') %>%
  unnest(mrd) %>%
  left_join(summary_tbl) %>%
  mutate(height_prior = mrd - TMRCA_lower) %>% 
  mutate(sg = ceiling(height_prior) * 6) %>%
  mutate(across(starts_with('TMRCA'), .fns = ~ subtract(mrd, .x), .names = 'height{.col}')) %>%
  view()


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################

summary_tbl %>%
  mutate(across(starts_with('TMRCA'), .fns = ~decimal_date(.x))) %>%
  mutate(clade = str_split_i(clade, '_', 3)) %>% 
  left_join(heights %>% select(clade, mrd)) %>%
  mutate(height_prior = mrd - TMRCA_lower) %>% 
  mutate(sg = ceiling(height_prior) * 6) %>% 
  mutate(across(starts_with('TMRCA'), 
                .fns = ~ subtract(mrd, .x), 
                .names = 'height{.col}')) %>% 
  view()