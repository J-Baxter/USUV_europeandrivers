################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-12-11
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


GetPerBranchMetrics <- function(localTreesDirectory){
  require(fields)
  
  tree_extractions <- list.files(path = localTreesDirectory,
                                 pattern = 'TreeExtractions_',
                                 full.names = T) %>%
    read_csv()
  
  out <- tree_extractions %>%
    bind_rows(.id = 'tree') %>%
    mutate(time = endYear - startYear,
           r1 = cbind(startLon, startLat),
           r2 = cbind(endLon, endLat)) %>%
    rowwise() %>%
    mutate(dists = fields::rdist.earth.vec(r1, r2, miles = F),
           bearing = geosphere::bearing(r1, r2),
           vel = (dists/time)/365,
           coef =( dists**2)/(4*time)) %>%
    ungroup()
  
  return(out)
}
################################### DATA #######################################
# Read and inspect data
per_branch_dirs <- dirs[grepl('Partial', dirs)]

per_branch_tbl <- paste(per_branch_dirs, '/Extracted_trees', sep = '') %>%
  lapply(., GetPerBranchMetrics) %>%
  setNames(per_branch_dirs) %>%
  bind_rows(.id = 'filename') %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade')) %>%
  dplyr::select(-starts_with('file'))


################################### MAIN #######################################
# Main analysis or transformation steps
per_branch_tbl %>%
  mutate(bearing = if_else(bearing < 0, bearing + 360, bearing))  %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>% 
  group_by(clade) %>%
  mutate(bearing = cut(x = bearing, 
                       breaks = seq(0, 360, by = 10),
                       labels = NULL)) %>%
  group_by(bearing, .add = T) %>%
  summarise(n = n(),
            wdc = mean(coef)) %>% 
  ggplot() +
  geom_col(aes(x = bearing, y= n, fill = wdc)) + 
  scale_y_continuous(NULL) +
  scale_x_discrete(NULL,
                   breaks = c("(0,10]", "(90,100]", "(180,190]", "(270,280]"),
                   labels = c("0", "90", "180", "270")) +
  scale_fill_distiller(trans = "log10", direction = 1) +
  coord_radial(rotate_angle = TRUE, expand = FALSE) + 
  facet_wrap(~clade, axes = 'all_x', scales = 'free_y', nrow = 3) + 
  theme_minimal(base_size = 14) + 
  theme(plot.margin=grid::unit(c(5,5,5,5), "mm"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold',  size = 14,margin = margin(1,0,1,0, "mm"))) +
  guides(fill= guide_colourbar('Weighted Diffusion Coefficient', 
                               position = 'inside',
                               theme = theme(legend.position.inside = c(0.8,0.15),
                                             legend.justification = c(0,0))))



################################### OUTPUT #####################################
# Save output files, plots, or results
ggsave('./2025Jun24/plots/phylogeo_direction.pdf',
       dpi = 360,
       height = 15,
       width = 10)
#################################### END #######################################
################################################################################