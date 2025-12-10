################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
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

ExtractRoots <- function(tree_list){
  out <- lapply(tree_list, function(x)  as_tibble(x) %>% 
                  filter(parent == node) %>%
                  pull(location) %>%
                  unlist() %>%
                  as.numeric()) %>%
    enframe() %>%
    unnest_wider(value, names_sep = '_') %>%
    rename(lat = value_1,
           lon = value_2) %>%
    st_as_sf(coords = c('lon', 'lat'),
             crs = st_crs(4326)) %>%
    st_transform(3035) %>%
    st_coordinates() %>%
    as_tibble()
  
  return(out)
}


################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]

treelog_files <- sapply(dirs, list.files, pattern = 'joint_500.trees|empirical_cont_500.trees',
                        full.names = TRUE, 
                        simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()

################################### MAIN #######################################
# Main analysis or transformation steps

treelogs <- lapply(treelog_files, read.beast) 

root_locations <- lapply(treelogs, ExtractRoots) %>%
  setNames(treelog_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'file4')) %>%
  separate_wider_delim(clade, delim = '_', names = c('data','clade')) %>%
  mutate(data = if_else(data == 'Partial', 'Combined', data)) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)")))



map <- ne_countries(continent = 'europe', scale = "medium", returnclass = "sf") %>%
  st_transform(3035)

ggplot() +
  geom_sf(data = map) +
  
  # Plot dots
  geom_hdr(data = root_locations, aes(x = X, y = Y, fill = factor(data, levels = c('Combined', 'NFLG'))),
           probs = c( 0.95, 0.75, 0.5)) +
  coord_sf(datum = sf::st_crs(3035),
           #ylim = c(1200000, 5500000),
           #xlim = c(2500000, 6700000), 
           ylim = c(2000000, 3500000),
           xlim = c(3500000, 5400000), 
           expand = FALSE) +
  scale_x_continuous(expand = c(0.01,0.01)) + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = c('Combined' = '#e41a1c', 'NFLG' = '#377eb8')) + 
  facet_wrap(~clade,
             nrow = 3) +
  theme_void(base_size = 14) +
  theme(plot.margin=grid::unit(c(5,5,5,5), "mm"),
        panel.spacing = unit(2.5, "lines"),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold',  size = 14,margin = margin(1,0,1,0, "mm"))) +
  guides(alpha= guide_legend('Highest Posterior Density', position = 'inside', theme = theme(legend.position.inside = c(0.75,0.15),
                                                                                             legend.justification = c(0,0))) ,
         fill= guide_legend('Data', position = 'inside', theme = theme(legend.position.inside = c(0.75,0.15),
                                                                       legend.justification = c(0,0))))


ggsave('./2025Jun24/plots/phylogeo_root_maps.pdf',
       dpi = 360,
       height = 15,
       width = 10)

################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################

