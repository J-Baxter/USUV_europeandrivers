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




FormatPhyloGeo <- function(mcc_file, posterior_file){
  require(rnaturalearthdata)
  require(tidytree)
  require(sf)
  require(tidyverse)
  require(seraphim)
  
  # Import tree data
  mcc_tree <- read.beast(mcc_file)
  tree_tbl <- as_tibble(mcc_tree)
  
  # Guess most_recent_date
  most_recent_date <- tree_tbl %>%
    slice_min(height_median, n = 1, with_ties = FALSE) %>%
    pull(label) %>%
    str_extract(., '\\d{4}-\\d{2}-\\d{2}') %>%
    ymd() %>%
    decimal_date() %>%
    unique()
  
  
  # Extract TMRCA
  start_date <- tree_tbl %>%
    slice_max(height_median, n = 1, with_ties = FALSE) %>% 
    pull(height_median) %>% 
    as.numeric() %>% 
    subtract(most_recent_date,.)
  print(start_date)
  
  # Scan posterior tree file
  allTrees <- scan(file = posterior_file,
                   what = '',
                   sep = '\n',
                   quiet = T)
  
  
  localTreesDirectory = "./2025Jun10/temp_tree_dir/"
  do.call(file.remove, list(list.files("./2025Jun10/temp_tree_dir/", full.names = TRUE)))
  
  burnIn <- 0
  randomSampling <- TRUE
  nberOfTreesToSample <- 100
  #mostRecentSamplingDatum <- most_recent_date
  coordinateAttributeName <- "location"
  treeExtractions(localTreesDirectory,
                  allTrees,
                  burnIn, 
                  randomSampling, 
                  nberOfTreesToSample, 
                  most_recent_date,
                  coordinateAttributeName,
                  nberOfCores = 8)
  
  # Step 4: Estimating the HPD region for each time slice ----
  
  # Format Polygons
  polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, 
                                              nberOfExtractionFiles = nberOfTreesToSample,
                                              prob = 0.95, 
                                              startDatum = start_date, 
                                              precision = 0.08))
  
  # Conver polgons to GEOMETRY, set coordinate system and add values
  polygons_sf <- lapply(polygons, st_as_sf)  %>%
    bind_rows() %>%
    pivot_longer(cols = starts_with('2'), names_to = 'year', values_to = 'value') %>%
    filter(value == 1) %>%
    mutate(year = as.numeric(year)) %>%
    dplyr::select(-value) %>%
    st_set_crs(4326) 
  
  
  # Format Nodes
  nodes_sf <- tree_tbl %>%
    dplyr::select(node,height_median, location1, location2, host_simplifiedhost, label) %>%
    mutate(height_median= as.numeric(height_median)) %>%
    replace_na(list(height_median = 0)) %>%
    mutate(year = most_recent_date - height_median) %>%
    
    # Convert to POINT & set coordinate system
    st_as_sf(coords = c( 'location2', 'location1'), 
             crs = 4326)
  
  # Format Arrows
  edges <-  tree_tbl %>%
    dplyr::select(node, parent) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2)) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2),
              by = join_by(parent== node)) %>%
    mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
    mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                   .default = location1.y),
           location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                   .default = location2.y)) %>%
    rename(start_lat = location1.x,
           start_lon = location2.x,
           end_lat = location1.y,
           end_lon = location2.y)
  
  edges_sf <- edges %>%
    rowid_to_column(var = 'id') %>%
    pivot_longer(starts_with(c('start', 'end')),
                 names_to = c("type", ".value"),
                 names_sep = "_") %>%
    
    # Convert coordinate data to sf POINT
    st_as_sf( coords = c("lon", "lat"),
              crs = 4326) %>%
    
    # Convert POINT geometry to MULTIPOINT, then LINESTRING
    group_by(id) %>% 
    summarise(do_union = FALSE) %>% 
    st_cast("LINESTRING") %>% 
    
    # Convert rhumb lines to great circles
    st_segmentize(units::set_units(20, km)) %>%
    
    # Wrap dateline correctly
    st_wrap_dateline()
  
  
  out <- list(polygons = polygons_sf,
              edges = edges_sf,
              nodes = nodes_sf)
  
  
  return(out)
}


PlotPhyloGeo <- function(phylogeo_list){
  polygons_sf <- phylogeo_list[['polygons']]
  edges_sf <- phylogeo_list[['edges']]
  nodes_sf <- phylogeo_list[['nodes']]
  
  # load base map
  map <- ne_countries(continent = 'europe', scale = "medium", returnclass = "sf")
  
  # Plot in GGplot
  plot <- ggplot() +
    geom_sf(data = map) +
    
    # Plot HPD polygons
    geom_sf(data = polygons_sf, 
            aes(fill = year), 
            lwd = 0, 
            alpha = 0.04) + 
    
    # Plot Branches
    geom_sf(data = edges_sf,
            lwd = 0.2) + 
    
    # Plot nodes
    geom_sf(data = nodes_sf,
            size = 1.5, 
            aes(fill = year, colour = year, shape = is.na(label)))+
    
    scale_shape_manual(values = c(19,1),
    ) + 
    
    # Set graphical scales - must be fixed
    scale_fill_viridis_c('Year',
                         limits = c(1998, 2025),
                         #breaks=c(2019, 2020,2021,2022,2023,2024),
                         direction = -1,
                         option = 'C')+
    
    scale_colour_viridis_c('Year',
                           limits = c(1998, 2025),
                           #breaks=c(2019, 2020,2021,2022,2023,2024),
                           direction = -1,
                           option = 'C')+
    
    coord_sf(ylim = c(-60, 80),
             xlim = c(-185, 185),
             expand = TRUE) +
    
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    
    guides(colour = guide_colourbar(
      theme = theme(
        legend.key.height  = unit(0.75, "lines"),
        legend.key.width = unit(10, "lines")),
      #title.position = 'left',
      title.vjust = 1,
      position = 'bottom'), 
      
      fill = guide_colourbar(
        theme = theme(
          legend.key.height  = unit(0.75, "lines"),
          legend.key.width = unit(10, "lines")),
        title.vjust =1,
        #title.position = 'left',
        position = 'bottom'), 
      
      shape = 'none') +
    
    theme_void() + 
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
          legend.text = element_text(size = 8),
          legend.position = 'none', 
          panel.spacing = unit(2, "lines"), 
          strip.background = element_blank()) 
  
  return(plot)
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



# by year
plt_c <- diffusioncoef_tbl %>%
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