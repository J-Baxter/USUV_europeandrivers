tree_logs %>%
  filter(grepl('^meanRate$', Parameter)) %>%
  rename_with(~gsub('^value', 'evolutionaryrate', .x) %>%
                gsub('\\.', '_', .),
              .cols = starts_with('value')) %>%
  mutate(clade = str_split_i(clade, '_', 3)) %>%
  mutate(include = ifelse(clade %in% c('IV', 'II'), 'Exclude', 'Include') %>% factor(levels = c('Include', 'Exclude'))) %>%
  ggplot(aes(y = clade, x = evolutionaryrate, slab_fill = clade)) + 
  stat_halfeye(p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95) +
  theme_minimal()+
  scale_y_discrete('Clade')+
  scale_x_continuous('Evolutionary Rate (Subs/Site/Year)',
                     breaks = seq(0.000, 0.0015, by = 0.00025)) +
  scale_fill_brewer(aesthetics = 'slab_fill', palette = 'Dark2') + 
  facet_grid(rows = vars(include), scales = 'free_y', space = 'free', switch = 'y')+
  theme(legend.position = 'none',
        strip.placement  = 'outside',
        strip.text = element_text(face = 'bold', size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))


ggsave('./2025Jun24/plots/withineurope_nflg_evolutionaryrates.jpeg',
       dpi = 360,
       height = 12,
       width = 17,
       units = 'cm')
  