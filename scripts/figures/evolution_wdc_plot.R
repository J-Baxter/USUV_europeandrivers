################################################################################
## Script Name:        Evolution of Weighted Diffusion Coefficient Plot
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-12-10
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


################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]

diffusioncoef_files <- sapply(dirs, list.files, pattern = 'tree_weighted_diffusion_coefficient.txt',
                              full.names = TRUE, 
                              simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()


################################### MAIN #######################################
# Main analysis or transformation steps
diffusioncoef_tbl <- lapply(diffusioncoef_files, ReadSpatialStatistics) %>%
  setNames(diffusioncoef_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(metric = gsub('USUV_tree_|\\.txt', '', metric))

diffusioncoef_tbl %>%
  mutate(draw = replace_na(draw,0)) %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  group_by(clade, time) %>%
  median_hdci(draw) %>%
  drop_na(draw) %>%
  ggplot() +
  geom_path(aes(x = time,
                y = draw), 
            linetype = 'dashed', 
            linewidth = 0.5) + 
  
  geom_ribbon(aes(x = time, 
                  ymin = .lower,
                  ymax = .upper),
              colour = '#ff000000',
              alpha = 0.45) + 
  
  facet_wrap(.~clade,
             ncol = 5,
             axes = 'all',
             scales = "free_y"
             # axes = 'all',
  )  +
  
  scale_y_continuous(expression(paste('Weighted Diffusion Coefficient (',Km^{2}, ' ', year^{-1}, ')' )), 
                     expand = c(0,0))+
  scale_x_continuous('Time', expand = c(0.01,0), 
                     limits = c(2000,2025),
                     breaks = seq(2000,2025, by = 5),
                     labels = seq(2000,2025, by = 5)) + 
  #coord_cartesian(ylim = c(0, 7500)) +
  
  theme_classic(base_size = 9) + 
  
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.8,
                                   hjust = 0.8),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold', size = 9)) 
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################