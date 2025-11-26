traits_tree <- read.beast('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits_mcc.tree')

traits_tree_tbl <-  read.beast('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits_mcc.tree') %>%
  as_tibble()
all_logs <- beastio::readLog('./2024Oct20/test_beast/traits/USUV_2024Oct20_nflg_subsample1_traits.log') %>%
  ggs()

shading_intervals <- seq(1920, 2030, by = 10)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2030)

p1 <- traits_tree %>% 
  ggtree(mrsd = most_recent_date, aes(colour = dist_50)) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20)) +
  scale_colour_brewer('Location', palette = 'Dark2', labels = c('not_europe' = 'Outside of Europe',
                                                                'europe' = 'Within Europe'))+
  new_scale_colour()+
  geom_nodepoint(aes(colour = dist_50))+
  scale_colour_brewer('Location', palette = 'Dark2', labels = c('not_europe' = 'Outside of Europe',
                                                                'europe' = 'Within Europe'),
                      guide="none") + 
  #Background
  annotate("rect", 
           xmin = shading_intervals[seq(1, length(shading_intervals), 2)], 
           xmax = shading_intervals_end[seq(1, length(shading_intervals_end), 2)], 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
  
  geom_density(data = tmrca, 
               aes(x = value, y = (..count../Ntip(all_europe))),
               bounds = GetNodeAges("2023-09-23",  c(70.5184216919521, 179.175463927162)) %>% decimal_date() %>% sort(),
               fill = "blue", 
               alpha = 0.5, 
               inherit.aes = FALSE)


p <- traits_tree %>% 

  ggtree(mrsd = most_recent_date) + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 20)) +
  
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
                     'New Sequences') +
  scale_colour_manual(values = c("TRUE" = 'blue'), 
                      'New Sequences') + 
  
  # node colour to show pp support
  new_scale_colour()+
  geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support') + 
  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = as.factor(is_europe)),
             width = 3,
             colour = "white",
             pwidth = 1,
             offset = 0.03) + 
  scale_fill_brewer('Location', palette = 'Dark2', labels = c('0' = 'Outside of Europe',
                                                              '1' = 'Within Europe')) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = nuts0_id),
             width = 3,
             colour = "white",
             pwidth = 1,
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
                           'RS' = 'Serbia',
                           'SE' = 'Sweden',
                           'SK' = 'Slovakia',
                           'SN' = 'Senegal',
                           'UG' = 'Uganda',
                           'UK' = 'United Kingdom',
                           'ZA' = 'South Africa' )) 