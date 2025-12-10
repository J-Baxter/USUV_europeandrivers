################################################################################
## Script Name:        Extract posterior estimates from NFLG trees
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
library(treeio)
library(beastio)
library(ggmcmc)
library(ggdist)

################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]


logfilepaths <- sapply(dirs, 
                       list.files,
                       pattern = "HMC_(SG|constant)\\.log$|thin.log|III_test.log",
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
  bind_rows(., .id = 'clade') %>% 
  separate_wider_delim(clade, delim = '/', names = c('file0', 'file1', 'file2', 'clade')) %>%
  separate_wider_delim(clade, delim = '_', names = c('data','clade')) %>%
  mutate(data = if_else(data == 'Partial', 'Combined', data)) %>%
  dplyr::select(-starts_with('file'))

################################### MAIN #######################################
# Main analysis or transformation steps
tree_logs %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  filter(Parameter == 'age.root.') %>%
  ggplot(aes(fill = data, x = value, colour = data)) + 
  stat_halfeye(point_interval = NULL,
               p_limits = c(0.01, 0.99),
               slab_alpha = 0.6) +
  stat_pointinterval(point_interval = "median_hdci",
                     y = rep(c(0.8,0.9), 5),
                     .width =  0.95) +
  scale_x_continuous('Time' , limits = c(1995, 2010), expand = c(0,0)) + 
  scale_y_continuous('Probability Density', expand = c(0.01,0), limits = c(0,1)) + 
  scale_fill_manual(values = c('Combined' = '#e41a1c', 'NFLG' = '#377eb8')) + 
  scale_colour_manual(values = c('Combined' = '#e41a1c', 'NFLG' = '#377eb8')) + 
  facet_wrap(~clade,
             axes = 'all',
             nrow = 3) +
  theme_classic(base_size = 14) + 
  theme(plot.margin=grid::unit(c(5,5,5,5), "mm"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold',  size = 14,margin = margin(1,0,1,0, "mm"))) +
  guides(fill= guide_legend('Data', position = 'inside', theme = theme(legend.position.inside = c(0.75,0.15),
                                                                       legend.justification = c(0,0))),
         colour= guide_legend('Data', position = 'inside', theme = theme(legend.position.inside = c(0.75,0.15),
                                                                       legend.justification = c(0,0))))

 


################################### OUTPUT #####################################
# Save output files, plots, or results
ggsave('./2025Jun24/plots/tmrca_all_data.pdf',
       dpi = 360,
       height = 15,
       width = 10)


#################################### END #######################################
################################################################################