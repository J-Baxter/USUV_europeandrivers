####################################################################################################
####################################################################################################
## Script name: NFLG MCC USUV tree visualisation
##
## Purpose of script:
##
## Date created: 2024-10-25
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 
options(ignore.negative.edge=TRUE)
  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(treeio)
library(TreeTools)
library(ggtree)
library(ggnewscale)
library(ggsci)
library(ggtreeExtra)

# User functions


############################################## DATA ################################################
nflg_mcc <- read.beast('./2025Feb10/global_analysis/beast_plain/USUV_2025Feb10_NFLG_SRD06_relaxLn_constant_mcc.tree')
nflg_ca <- read.beast('./2025Feb10/global_analysis/beast_plain/USUV_2025Feb10_NFLG_SRD06_relaxLn_constant_ca.tree')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2025Feb10.csv') %>%
  filter(tipnames %in% nflg_ca@phylo$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(nflg_mcc@phylo)) #sanity check


############################################## MAIN ################################################

most_recent_date <- metadata_in_tree %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym
  max(na.rm = TRUE)

#ggtree(nflg_ca,
      # mrsd = most_recent_date) + 
 # theme_tree2() +

p <- nflg_ca %>% 
  
  left_join(metadata_in_tree %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames, lineage, host_class) %>%
              rename(label = tipnames),
            by = 'label') %>%
  ggtree(mrsd = most_recent_date) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 16)) +
  
  scale_x_continuous(
    #limits = c(2000, 2023),
    breaks = seq(1920, 2020, 20)) +
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = is.na(sequence_accession),
                    shape = is.na(sequence_accession)),
                size = 3, 
                alpha = 0.9) +
  
  geom_tiplab(aes(colour = is.na(sequence_accession)),
              #align = TRUE, 
              size = 0) +

  scale_shape_manual(values = c("TRUE" = 18),
                     'New Sequences',
                     guide = 'none') +
  scale_colour_manual(values = c("TRUE" = 'blue'), 
                      'New Sequences',
                      guide = 'none') + 

  
  # node colour to show pp support
  new_scale_colour()+
  geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support',
                        guide = guide_colourbar(order = 4)) + 
  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = as.factor(is_europe)),
             width = 4,
             colour = "white",
             pwidth = 1.2,
             offset = 0.03) + 
  scale_fill_brewer('Location', palette = 'Dark2', direction = -1, 
                    labels = c('0' = 'Outside of Europe',
                               '1' = 'Within Europe'),
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 2)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = nuts0_id),
             width = 4,
             colour = "white",
             pwidth = 1.2,
             offset = 00.04) +
  scale_fill_d3(name = 'Country', 
                palette ='category20', 
                alpha = 0.99, 
                labels = c('AT' = 'Austria',
                           'BE' = 'Belgium',
                           'CF' = 'Central African Republic',
                           'DE' = 'Germany',
                           'ES' = 'Spain',
                           'FR' = 'France',
                           'HU' = 'Hungary',
                           'IL' = 'Israel',
                           'IT' = 'Italy',
                           'LU' = 'Luxembourg',
                           'NL' = 'Netherlands',
                           'PT' = 'Portugal',
                           'RS' = 'Serbia',
                           'SE' = 'Sweden',
                           'SK' = 'Slovakia',
                           'SN' = 'Senegal',
                           'UG' = 'Uganda',
                           'UK' = 'United Kingdom',
                           'EL' = 'Greece',
                           'ZA' = 'South Africa' ),
                guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 2, order = 1)) +
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = host_class),
             width = 4,
             colour = "white",
             pwidth = 1.2,
             offset = 00.04) +
  scale_fill_brewer('Host Class', 
                    labels = c('Aves' = 'Bird',
                               'insecta'= 'Insect',
                               'mammalia' = 'Mammal',
                               'arachnida' = 'Arachnid'),
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 3))+
  theme(legend.position = c(0.2,0.6),
       # legend.position = "bottom",       # Place legends at the bottom
        legend.box = "vertical",
       legend.direction = 'vertical'
        )




ggsave(p, filename = "~/Downloads/summarytree.png", width = 13, height = 12)

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################