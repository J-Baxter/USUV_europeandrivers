################################################################################
## Script Name:        Figure 2 Composite
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-11-27
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
## Note this script requires sections from nflg_patristicdistances.R and
# cluster_evaluation.R

CladeLab <- function(fort_tree, scale = NULL){
  # Determine vertical extent of clade
  
  ymin <- min(fort_tree$y[fort_tree$isTip])
  ymax <- max(fort_tree$y[fort_tree$isTip])
  dif <- ymax - ymin
  
  if(is.numeric(scale)){
    ymin <- ymin-scale*dif
    ymax <- ymax+scale*dif
  }
  
  # Determine root x position (or any x for the bar)
  x_bar <- max(fort_tree$x) + 1   # offset slightly left
  
  out <- tibble('ymin' = ymin,
                'ymax' = ymax,
                'x_bar'  = x_bar)
  return(out)
  
}


################################### DATA #######################################
# Read and inspect data
nflg_hipstr <- read.beast('./2025Jun24/global_analysis/USUV_2025Jun24_NFLG_SRD06_HMC_empirical_hipstr.tree')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')%>%
  filter(tipnames %in% nflg_hipstr@phylo$tip.label) 

metadata <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')
################################### MAIN #######################################
# Main analysis or transformation steps

a <-  ggtree(nflg_hipstr, 
               mrsd = '2024-10-12',
               aes(colour = is_europe)) + 
  scale_colour_manual( 'Continent' ,
                       values = c( '#2f91bd', '#24a489'), 
                     labels = c('Europe', 'Non-Europe')) + 
  theme_tree2(base_size = 10) +
  scale_y_reverse(expand = expansion(mult = c(0, .05))) + 
  theme(legend.position = 'inside',
        legend.background = element_blank(),
        legend.position.inside = c(0,0.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.justification=c(0,1),
        text = element_text(size = 9))


b <- metadata %>%
  count(is_europe) %>%
  mutate(is_europe = if_else(is_europe == 1, 'Europe', 'Non-Europe')) %>%
  ggplot(aes(x = is_europe, fill = is_europe, y = n, colour = is_europe)) + 
  geom_bar(stat = 'identity', alpha = 0.9) + 
  scale_fill_manual(values = c( '#2f91bd', '#24a489'), 
                    labels = c('Europe', 'Non-Europe')) + 
  scale_colour_manual(values = c( '#2f91bd', '#24a489'), 
                      labels = c('Europe', 'Non-Europe')) + 
  theme_classic(base_size =9)+
  scale_x_discrete('Continent', expand = c(0.3,0.3)) + 
  scale_y_continuous('N', expand = c(0,0)) + 
  theme(legend.position = 'null')


c <-  bind_rows(within_europe,
                outside_europe_within,
                outside_europe_between) %>%
  ggplot() +
  geom_density(aes(x = distance, fill = type, colour = type),alpha = 0.5)+
  scale_fill_manual(NULL, values = c(  '#24a489', 'darkgrey','#2f91bd')) + 
  scale_colour_manual(NULL, values = c( '#24a489', 'darkgrey',  '#2f91bd')) + 
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic(base_size = 9) +
  guides(fill=guide_legend(nrow=3),
         colour=guide_legend(nrow=3))


# import clusters
all_clusters <- read_csv('./2025Jun24/europe_clusters/all_clusters.csv') %>%
  group_by(cluster) %>%
  dplyr::mutate(cluster = dplyr::cur_group_id()) %>%
  ungroup() %>%
  mutate(cluster = as.roman(cluster))

# comman ancestor of tips in tree
common_ancestors <- nflg_hipstr %>%
  as_tibble() %>%
  left_join(all_clusters) %>%
  split(~ as.character(cluster)) %>%
  lapply(., pull, node) %>%
  lapply(., MRCA, .data = nflg_hipstr) %>%
  flatten_dbl() 

lineage_nodes <- common_ancestors %>%
  offspring(nflg_hipstr, ., self_include = F) %>% 
  flatten_dbl() 

# Extract subtrees
sub1 <- tree_subset(nflg_hipstr, 503, levels_back = 0, group_name = 'sub1')
sub2 <-tree_subset(nflg_hipstr, 504, levels_back = 0, group_name = 'sub2')
sub3 <-tree_subset(nflg_hipstr, 514, levels_back = 0, group_name = 'sub3')
sub4 <-tree_subset(nflg_hipstr, 710, levels_back = 0, group_name = 'sub4')
sub5 <-tree_subset(nflg_hipstr, 713, levels_back = 0, group_name = 'sub5')
sub6 <-tree_subset(nflg_hipstr, 834, levels_back = 0, group_name = 'sub6')
sub7 <-tree_subset(nflg_hipstr, 815, levels_back = 0, group_name = 'sub7')
sub8 <-tree_subset(nflg_hipstr, 811, levels_back = 0, group_name = 'sub8')

# Convert to ggtree-compatible data frames
d1 <- fortify(sub1) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) 

d2 <- fortify(sub2) %>% 
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d1$y) + 5) 

d3 <- fortify(sub3) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d2$y) + 7) 

d4 <- fortify(sub4) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d3$y) + 7)

d5 <- fortify(sub5) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d4$y) + 15)

d8 <- fortify(sub8) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d5$y) +15)

d7 <- fortify(sub7) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d8$y) + 22)

d6 <- fortify(sub6) %>%
  mutate(x = (2024.779 -as.numeric(height_median))) %>%
  mutate(y = y + max(d7$y)+ 22)


# Base plot: empty coordinate system
d <- ggplot() +
  geom_tree(data = d1, aes(x = x,
                           y = y),
            color =  '#2f91bd') + #1
  geom_tree(data = d2, aes(x = x,
                           y = y),
            color =  '#2f91bd') + #2
  geom_tree(data = d3, aes(x =x, 
                           y = y), 
            color =  '#2f91bd') + #3
  geom_tree(data = d4, aes(x = x,
                           y = y), 
            color =  '#2f91bd') + #4
  geom_tree(data = d5, aes(x = x,
                           y = y), 
            color =  '#2f91bd') + #5
  geom_tree(data = d8, aes(x = x,
                           y = y),
            color =  '#2f91bd') + #8
  geom_tree(data = d7, aes(x = x,
                           y = y),
            color =  '#2f91bd') + #7
  geom_tree(data = d6, aes(x = x,
                           y = y), 
            color =  '#2f91bd') + #6
  
  theme_tree2()  +
  scale_y_reverse(expand = expansion(mult = c(0.01, .01))) + 
  new_scale_colour()+
  geom_tippoint(data = d1[d1$isTip, ], 
                aes(x = x, y = y ), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d2[d2$isTip, ], aes(x = x,y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d3[d3$isTip, ], aes(x = x, y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d4[d4$isTip, ], aes(x = x,y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d5[d5$isTip, ], aes(x = x,y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d8[d8$isTip, ], aes(x = x, y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d7[d7$isTip, ], aes(x = x, y = y), 
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_tippoint(data = d6[d6$isTip, ], aes(x = x, y = y),
                colour = 'black',
                fill = '#2f91bd',
                size = 2, 
                shape = 21) +
  
  geom_segment(data = CladeLab(d1, scale = 2.5),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d1),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "I",
            size = 3) + 
  
  geom_segment(data = CladeLab(d2, scale = 1.5),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d2),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "II",
            size = 3) + 
  
  geom_segment(data = CladeLab(d3),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d3),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "III",
            size = 3) +
  
  geom_segment(data = CladeLab(d4, scale = 2),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d4),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "IV",
            size = 3) + 
  
  geom_segment(data = CladeLab(d5),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d5),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "V",
            size = 3) + 
  
  geom_segment(data = CladeLab(d6),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d6),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "VI",
            size = 3) +
  
  geom_segment(data = CladeLab(d7),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d7),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "VII",
            size = 3)+ 
  
  geom_segment(data = CladeLab(d8, scale = 2),
               aes(x = x_bar, xend = x_bar, y = ymin, yend = ymax),
               color = "black", size = 0.7)  +
  geom_text(data = CladeLab(d8),
            aes(x = x_bar + 1, 
                y = (ymin + ymax)/2),
            label = "VIII",
            size = 3) +
  theme(text = element_text(size = 9))



# Elbow curve
e <- cluster_verylong %>%
  summarise(n = n_distinct(cluster), .by = distance_threshold) %>%
  ggplot() + 
  geom_vline(aes(xintercept = 35), linetype = 'dashed', colour = '#2f91bd')+
  geom_point(aes(y = n, x = as.numeric(distance_threshold))) +
  scale_y_continuous('Clusters (n)') + 
  scale_x_continuous('Patristic Distance Threshold', expand = c(0,0)) +
  theme_classic(base_size = 9) 


# Patristic distance plot
f <- tmp_2 %>%
  pivot_longer(c(starts_with('type') ),
               names_to = 'cluster_threshold',
               values_to = 'cluster_type') %>%
  mutate(cluster_threshold = gsub('type_', '', cluster_threshold) %>%
           as.numeric()) %>%
  select(pair_id, distance, label, starts_with('cluster')) %>%
  filter(cluster_threshold == 35) %>%
  ggplot() +
  geom_histogram(aes(x = distance, 
                     y=..density..,
                     fill = cluster_type, 
                     colour = cluster_type),
                 alpha = 0.5,
                 binwidth = 2.5,
                 position="identity" )+
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL) + 
  scale_fill_d3(palette = 'category20',
                alpha = 0.99,
                name= NULL) +
  scale_x_continuous('Patristic Distance',
                     expand = c(0,0))+
  scale_y_continuous('Probability Density',
                     expand = c(0.00,0))+
  theme_classic(base_size = 9) +
  guides(fill=guide_legend(nrow=2, position = 'top'),
         colour=guide_legend(nrow=2, position = 'top'))





################################### OUTPUT #####################################
# Save output files, plots, or results
c_legend <- get_legend(c) %>% as_ggplot() + theme(legend.margin=margin(c(0,0,0,0)))
rhs <- align_plots(c + theme(legend.position = 'none'), c_legend , e,f,d, align = 'hv', axis = 'rb' )
top <- plot_grid(a, b,rhs[[1]],  labels = 'AUTO', align = 'h', axis = 'tb', nrow = 1, label_size = 10, scale = 0.95)
right <- plot_grid(rhs[[2]], rhs[[3]], rhs[[4]], labels = c('', 'E', 'F'), align = 'v', axis = 'lr', ncol = 1, rel_heights = c(0.2,0.45,0.55), label_size = 10)

bottom <- plot_grid(rhs[[5]], right, labels = c('D', '', ''), ncol = 2, align = 'v', axis = 'tb', rel_widths = c(0.65,0.35), label_size = 10)
plot_grid(top, bottom, nrow = 2, align = 'v', axis = 'tblr', rel_heights = c(0.25,0.75))

ggsave('~/Downloads/figure2_mk2.pdf', dpi = 360, height = 25, width = 20, units = 'cm')
#################################### END #######################################
################################################################################