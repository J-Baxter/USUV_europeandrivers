################################################################################
## Script Name:       Plot spatial statistics and phylogeography 
## Purpose:           Plot phylogeography on a map, alongside spatial statistics
##                    including wavefront distances, diffusion coefficients and 
##.                   isolation-by-distance (IBD) signal
## Author:             James Baxter
## Date Created:       2025-11-12
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
# Load required libraries
library(tidyverse)
library(magrittr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(diagram)
library(lubridate)
library(treeio)
library(scales)
library(tidybayes)
library(ggdensity)
library(ggdist)

ReadSpatialStatistics <- function(filepath){
  out <-  read_delim(filepath,
                     delim = '\t', 
                     name_repair = function(nm) {
                       nm[2:length(nm)] <- paste0(seq_len(length(nm) - 1))
                       nm}) %>%
    pivot_longer(-1, names_to = 'iter', values_to = 'draw') 
  
  return(out)
  
}


GetPerBranchMetrics <- function(localTreesDirectory){
  require(fields)
  
  tree_extractions <- list.files(path = localTreesDirectory,
                                 pattern = 'TreeExtractions_',
                                 full.names = T) %>%
    read_csv()
  
  out <- tree_extractions %>%
    bind_rows(.id = 'tree') %>%
    mutate(time = endYear - startYear,
           r1 = cbind(startLat, startLon),
           r2 = cbind(endLat, endLon)) %>%
    rowwise() %>%
    mutate(dists = fields::rdist.earth.vec(r1, r2, miles = F),
           vel = (dists/time)/365,
           coef =( dists**2)/(4*time)) %>%
    ungroup()
  
  return(out)
}



################################### DATA #######################################
# Read and inspect data
dirs <- list.dirs("./2025Jun24/europe_clusters", recursive = FALSE)
dirs <- dirs[grepl("_III|_V|_VI", basename(dirs))]

wavefront_files <- sapply(dirs, list.files, pattern = 'wavefront_distances.txt',
                          full.names = TRUE, 
                          simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()

diffusioncoef_files <- sapply(dirs, list.files, pattern = 'tree_weighted_diffusion_coefficient.txt',
                          full.names = TRUE, 
                          simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()

summary_files <- sapply(dirs, list.files, pattern = 'estimated',
                              full.names = TRUE, 
                              simplify = F) %>%
  Filter(length,.) %>% 
  flatten_chr()


wavefront_tbl <- lapply(wavefront_files, ReadSpatialStatistics) %>%
  setNames(wavefront_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(metric = gsub('USUV_|s.txt', '', metric))

diffusioncoef_tbl <- lapply(diffusioncoef_files, ReadSpatialStatistics) %>%
  setNames(diffusioncoef_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-starts_with('file')) %>%
  mutate(metric = gsub('USUV_tree_|\\.txt', '', metric))
  
summary_tbl <- lapply(summary_files, read_delim) %>%
  setNames(summary_files) %>%
  bind_rows(.id = 'filename')  %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade', 'metric')) %>%
  dplyr::select(-c(starts_with('file'), metric))

per_branch_dirs <- dirs[grepl('Partial', dirs)]

per_branch_tbl <- paste(per_branch_dirs, '/Extracted_trees', sep = '') %>%
  lapply(., GetPerBranchMetrics) %>%
  setNames(per_branch_dirs) %>%
  bind_rows(.id = 'filename') %>%
  separate_wider_delim(filename, delim = '/', names = c('file0', 'file1', 'file2', 'clade')) %>%
  dplyr::select(-starts_with('file'))


################################### MAIN #######################################
# Main analysis or transformation steps
# Plots




# 1. Spatial Wavefront Distance
# 2. Patristic Wavefront Distance
plt_b <-wavefront_tbl %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  group_by(clade, metric, time) %>%
  median_hdci(draw) %>%
  ggplot() +
  geom_path(aes(x = time,
                y = draw,
                colour = metric), 
            linetype = 'dashed', 
            linewidth = 0.5) + 
  
  geom_ribbon(aes(x = time, 
                  ymin = .lower,
                  ymax = .upper, 
                  fill = metric),
              colour = '#ff000000',
              alpha = 0.45)+ 
  
  facet_grid(cols = vars(clade),
             axes = 'all')  +
  
  scale_y_continuous('Wavefront Distance (Km)', expand = c(0,0))+
  scale_x_continuous('Time', expand = c(0.01,0)) + 
  
  scale_color_discrete(labels = c('Patristic Wavefront Distance', 
                                  'Spatial Wavefront Distance')) + 
  
  scale_fill_discrete(labels = c('Patristic Wavefront Distance', 
                                 'Spatial Wavefront Distance')) +
  
  coord_cartesian(xlim = c(1998, 2025)) + 
  
  theme_classic(base_size = 9) + 
  
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.8,
                                   hjust = 0.8),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold', size = 9)) + 
  
  guides(fill= guide_legend(NULL, position = 'bottom'),
          colour = guide_legend(NULL, position = 'bottom'))


# 3. Weighted Diffusion Coefficient
# summary 
summary_tbl %>%
  ggplot(aes(x = weighted_diffusion_coefficient)) +
  stat_halfeye(p_limits = c(0.00001, 0.9999999),
               point_interval = "median_hdci",
               .width =  0.95
  ) +
  scale_x_continuous(limits = c(0,5500)) + 
  facet_grid(cols = vars(clade))


# 4. Diffusion Coefficient Variation 
plt_d <- summary_tbl %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  ggplot(aes(x = weighted_diffusion_coefficient, y = diffusion_coefficient_variation_among_branches)) + 
  geom_hdr(aes(fill = after_stat(probs))) + 
  facet_grid(~clade, axes = 'all') +
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km^{2}, ' ', year^{-1}, ')' )), 
                     expand = c(0,0))+
  scale_y_continuous('Diffusion coefficient variation among branches', expand = c(0.01,0)) +
  scale_fill_brewer() + 
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.8,
                                   hjust = 0.8),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold', size = 9)) +
  guides(fill= guide_legend('Highest Density Regions', position = 'bottom'),
         alpha = 'none')

stat_hdr_raster
scientific_10 <- function(x) {   parse(text=gsub("e+*", " %*% 10^", scales::scientific_format()(x)))}

plt_d <- per_branch_tbl %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  ggplot(aes(x = dists, y = coef)) +
  geom_hdr(aes(fill = after_stat(probs))) + 
  scale_fill_brewer() +
  scale_x_log10(expression(paste("Great Circle Distance", ' (', Log[10], ' Km)')),
                breaks = trans_breaks("log10", function(x) 10^x)(10**seq(-3, 3, by = 1.5)),
                labels = trans_format("log10", math_format(.x)),
                expand = c(0,0)) + 
  scale_y_log10(expression(paste("Branch Diffusion Coefficient", ' (', Log[10], ' ', Km**2, ' ',year**-1, ')')),
                breaks = trans_breaks("log10", function(x) 10^x)(10**seq(-5, 5, by = 2.5)),
                labels = trans_format("log10", math_format(.x)),
                expand = c(0,0))+
  theme_classic(base_size = 9) + 
  facet_grid(~clade, axes = 'all') +
  guides(fill= guide_legend('Highest Density Regions', position = 'bottom'),
         alpha = 'none') +
  theme(strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(face = 'bold', size = 9)) 




plot_grid(plt_b, plt_c, plt_d, nrow = 3, ncol = 1,  align = 'hv', axis = 'tblr', labels = 'AUTO', label_size = 10, scale = 0.95)


# 5. Isolate by distance
summary_tbl %>%
  mutate(clade = factor(clade, labels = c( "III (B)",
                                           "V (A.1/A.3)",
                                           "VI (A)",
                                           "VII (A.2)",
                                           "VIII (A/A.0)"))) %>%
  group_by(clade) %>%
  median_hdci(isolation_by_distance_signal_rP2)


summary_tbl %>%
  ggplot(aes(x = isolation_by_distance_signal_rP2, y= clade)) +
  stat_halfeye(p_limits = c(0.0001, 0.9999),
               point_interval = "median_hdci",
               .width =  0.95
  ) +
  facet_grid(rows = vars(clade), scales = 'free_y') + 
  scale_x_continuous(limits = c(0,1))

summary_tbl %>%
  ggplot(aes(x = isolation_by_distance_signal_rP1)) +
  stat_slab(density = ggdist::density_bounded(bandwidth = "nrd0", bounds = c(0, 1)), color = "gray25") +
  
  facet_grid(rows = vars(clade), scales = 'free_y')

################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################