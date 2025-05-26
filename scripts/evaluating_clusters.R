cluster_verylong <- lapply(seq(5, 80, by = 5),
                       InferClusters, 
                       phylo = nflg_ca@phylo,
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
  mutate(type_5 = case_when(n_distinct(dist_5) == 1 ~ 'Within Europe: Within Clade',
                          .default = 'Within Europe: Between Clade')) %>%
  mutate(type_10 = case_when(n_distinct(dist_10) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_51 = case_when(n_distinct(dist_15) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_20 = case_when(n_distinct(dist_20) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_25 = case_when(n_distinct(dist_25) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_30 = case_when(n_distinct(dist_30) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_35 = case_when(n_distinct(dist_35) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_40 = case_when(n_distinct(dist_40) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_45 = case_when(n_distinct(dist_45) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_50 = case_when(n_distinct(dist_50) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_55 = case_when(n_distinct(dist_55) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_60 = case_when(n_distinct(dist_60) == 1 ~ 'Within Europe: Within Clade',
                            .default = 'Within Europe: Between Clade')) %>%
  mutate(type_70 = case_when(n_distinct(dist_70) == 1 ~ 'Within Europe: Within Clade',
                             .default = 'Within Europe: Between Clade')) %>%
  mutate(type_80 = case_when(n_distinct(dist_80) == 1 ~ 'Within Europe: Within Clade',
                             .default = 'Within Europe: Between Clade'))


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
                     fill = type_80, colour = type_80),
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


