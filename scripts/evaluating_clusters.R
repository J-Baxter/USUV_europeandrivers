cluster_verylong <- lapply(seq(5, 80, by = 5),
                       InferClusters, 
                       phylo = beast_mcc@phylo,
                       n_threshold = 2, 
                       filter = TRUE,
                       metadata) %>%
  setNames(as.character(seq(5, 80, by = 5))) %>%
  bind_rows(., .id = 'distance_threshold') 


# Elbow method
cluster_verylong %>%
  summarise(n = n_distinct(cluster), .by = distance_threshold) %>%
  ggplot() + 
  geom_point(aes(y = n, x = as.numeric(distance_threshold)))


cluster_wide <- cluster_verylong %>%
  pivot_wider(values_from = cluster,
              names_from = distance_threshold,
              names_prefix = 'dist_')

most_recent_date <- '2024-10-12'

nflg_mcc %>% 
  left_join(cluster_wide) %>%
  ggtree(mrsd = most_recent_date) + 
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = dist_35))



pdist_clades <- lapply()
  
  
  test <- function(x, dist_mat){
    
  }
  
  
tmp <- patristic_distances %>%
  as.matrix() %>%
  .[europe_tips,europe_tips] %>%
  as_tibble(rownames = 'a') %>%
  pivot_longer(-a, names_to = 'b', values_to =  'distance') %>%
  filter(a != b) %>%
  rowid_to_column(var = 'pair_id') %>%
  pivot_longer(-c(distance, pair_id), values_to = 'label') %>%
    
  left_join(cluster_wide) 


tmp_2 <- tmp %>%
  group_by(pair_id) %>%
  mutate(across(
    starts_with("dist_"),
    ~ case_when(
      n_distinct(.x) == 1 ~ 'Within Europe: Within Clade',
      TRUE ~ 'Within Europe: Between Clade'
    ),
    .names = "type_{str_extract(.col, '\\\\d+')}"
  ))


tmp_3 <- tmp_2 %>%
  ungroup() %>%
  pivot_wider(id_cols = c(pair_id, distance, starts_with('type')), names_from = 'name', values_from =  c(label, starts_with('dist_'))) %>%
  filter(!duplicated(paste0(pmax(label_a, label_b), pmin(label_a, label_b))))
  


dists <- tmp_3 %>%
  dplyr::select(c(starts_with('type'), distance)) %>%
  mutate(across(starts_with('type'), .fns = ~  gsub('within europe\\: ', '', str_to_lower(.x)))) %>%
  pivot_longer(cols = starts_with('type'), names_to = 'threshold', values_to = 'group') %>%
  drop_na(group, threshold) %>%
  group_by(threshold, group) %>%
  reframe(stats = glance(fitdistr(distance, 'normal')))
  

tmp_data <- tmp_3 %>%
  dplyr::select(c(starts_with('type'), distance)) %>%
  mutate(across(starts_with('type'), .fns = ~  gsub('within europe\\: ', '', str_to_lower(.x)))) %>%
  pivot_longer(cols = starts_with('type'), names_to = 'threshold', values_to = 'group') %>%
  filter(threshold == 'type_35') %>%
  filter(group == 'between clade') %>%
  pull(distance)

fit <- fitdistr(tmp_data, 'normal')
glance(fit)
ggplot(tmp_3) +
  geom_histogram(aes(x = distance, 
                     #y = after_stat(ndensity),
                     fill = type_35, colour = type_35),
                 alpha = 0.5,
                 binwidth = 5,
                 position="identity" )+
  scale_colour_d3(palette = 'category20',
                  alpha = 0.99,
                  name= NULL) + 
  scale_fill_d3(palette = 'category20',
                alpha = 0.99,
                name= NULL) +
  scale_x_continuous('Patristic Distance',
                     expand = c(0.005,0),
                     limits = c(0,300))+
  scale_y_continuous('Probability Density',
                     expand = c(0.005,0))+
  theme_classic()+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.75))


