# sub trees


a_tree <-  read.beast('./2024Oct20/alignments/subset_alignments/run/USUV_2024Oct20_noFLI_NFLG_A_subsampled_mcc.tree')
b_tree  <-  read.beast('./2024Oct20/alignments/subset_alignments/run/USUV_2024Oct20_NFLG_B_subsampled_mcc.tree')
c_tree <-  read.beast('./2024Oct20/alignments/subset_alignments/run/USUV_2024Oct20_noFLI_NFLG_C_subsampled_mcc.tree')

a_most_recent_date <- metadata %>%
  left_join(patristic_distance_clusters,
            join_by(tipnames ==label)) %>%
  drop_na(dist_50) %>%
  filter(dist_50 == 'A') %>% 
  pull(date_ymd) %>%
  max(., na.rm = T)


b_most_recent_date <- metadata %>%
  left_join(patristic_distance_clusters,
            join_by(tipnames ==label)) %>%
  drop_na(dist_50) %>%
  filter(dist_50 == 'B') %>% 
  pull(date_ymd) %>%
  max(., na.rm = T)


c_most_recent_date <- metadata %>%
  left_join(patristic_distance_clusters,
            join_by(tipnames ==label)) %>%
  drop_na(dist_50) %>%
  filter(dist_50 == 'C') %>% 
  pull(date_ymd) %>%
  max(., na.rm = T)


a_tree %>% 
  
  left_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames, lineage) %>%
              rename(label = tipnames),
            by = 'label') %>%
  ggtree(mrsd = a_most_recent_date) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20)) +

  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = is.na(sequence_accession),
                    shape = is.na(sequence_accession)),
                size = 3, 
                alpha = 0.9) +
  
  geom_tiplab(aes(colour = is.na(sequence_accession)),
              #align = TRUE, 
              size = 0) +
  
  scale_shape_manual(values = c("TRUE" = 18),
                     'New Sequences') +
  scale_colour_manual(values = c("TRUE" = 'blue'), 
                      'New Sequences') + 
  
  # node colour to show pp support
  new_scale_colour()+
  geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support') + 

  geom_fruit(geom = geom_tile,
             mapping = aes(fill = nuts0_id),
             width = 1,
             colour = "white",
             pwidth = 1,
             offset = 00.06) +
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
                           'RS' = 'Serbia',
                           'SE' = 'Sweden',
                           'SK' = 'Slovakia',
                           'SN' = 'Senegal',
                           'UG' = 'Uganda',
                           'UK' = 'United Kingdom',
                           'ZA' = 'South Africa',
                           'CZ' = 'Czech Republic',
                           'RO' = 'Romania')) 


c_tree %>% 
  
  left_join(metadata %>% 
              dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames, lineage) %>%
              rename(label = tipnames),
            by = 'label') %>%
  ggtree(mrsd = c_most_recent_date) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20)) +
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = is.na(sequence_accession),
                    shape = is.na(sequence_accession)),
                size = 3, 
                alpha = 0.9) +
  
  geom_tiplab(aes(colour = is.na(sequence_accession)),
              #align = TRUE, 
              size = 0) +
  
  scale_shape_manual(values = c("TRUE" = 18),
                     'New Sequences') +
  scale_colour_manual(values = c("TRUE" = 'blue'), 
                      'New Sequences') + 
  
  # node colour to show pp support
  new_scale_colour()+
  geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support') + 
  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = nuts0_id),
             width = 1,
             colour = "white",
             pwidth = 1,
             offset = 00.06) +
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
                           'RS' = 'Serbia',
                           'SE' = 'Sweden',
                           'SK' = 'Slovakia',
                           'SN' = 'Senegal',
                           'UG' = 'Uganda',
                           'UK' = 'United Kingdom',
                           'ZA' = 'South Africa',
                           'CZ' = 'Czech Republic',
                           'RO' = 'Romania')) 

