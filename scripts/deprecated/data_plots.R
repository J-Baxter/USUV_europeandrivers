#data plots


# dependencies
library(terra)
library(sf)
library(rnaturalearth)
library(giscoR)
library(tidyverse)


# load data
metadata_noconcat <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')




# Plt - sequence length
b <- metadata_noconcat %>%
  ggplot()+
  #geom_rect(aes(xmin = 9000, xmax = Inf, ymin = 0, ymax = Inf), alpha = 0.007, fill = 'lightgrey')+
  geom_histogram(aes(x = sequence_length), binwidth = 150) +
  #geom_vline(xintercept  = 9000, linetype = 'dashed', linewidth = 1) + 
  geom_vline(xintercept  = 200, linetype = 'dashed', linewidth = 1) + 
  scale_x_continuous('Sequence Length', expand = c(0,0), breaks = seq(0,10000, by= 2000)) + 
  scale_y_continuous('Count', expand = c(0,0)) +
  scale_fill_brewer(palette = 'PuBu', na.value = 'black', 'Genomic Region') + 
  theme_classic(base_size = 9)


# Plt - alignment
c <- metadata_noconcat %>%
  filter(generegion_nflg == 0) %>%
  arrange(sequence_start) %>%
  dplyr::select(tipnames, sequence_start, sequence_end, nuts0_id) %>%
  pivot_longer(starts_with('sequence'), values_to = 'position', names_to = 'name') %>%
  dplyr::select(-name) %>%
  group_by(tipnames) %>%
  complete(position = seq(min(position), max(position))) %>%
  fill(nuts0_id, .direction = 'down') %>%
  fill(nuts0_id, .direction = 'down') %>%
  ungroup() %>%
  summarise(prop_coverage = n()/856, .by = position) %>%
  arrange(position) %>%
  ggplot() + 
  geom_ribbon(aes(x = position, ymax = prop_coverage, ymin = 0),
              colour = 'black',
              fill = 'blue',
              alpha = 0.5) +
  #geom_path(aes(x = position, y = prop_coverage), linewidth = 1) +
  scale_y_continuous('Proportion Coverage', expand = c(0,0)) + 
  scale_x_continuous('Alignment Position', expand= c(0.01, 0)) + 
  coord_cartesian(xlim = c(0,10000)) +
  theme_classic(base_size = 9)

#metadata_noconcat %>%
#  arrange(sequence_start) %>%
 # filter(generegion_nflg == 0) %>%
 # rowid_to_column(var = 'index') %>%
 # ggplot() + 
 # geom_linerange(aes(y = index, xmin = sequence_start , xmax = sequence_end), linewidth = 0.5) + 
 # scale_y_continuous(expand = c(0.01,0.01)) + 
 # scale_x_continuous(expand = c(0.05,0.01), name = 'Genome Coordinates', limits = c(0, 10500)) +
 # theme_classic(base_family = "LM Sans 10",
               # base_size = 20)


nongisco_shapefile <- ne_countries(scale = 10, country = c('bosnia and herzegovina', 'kosovo', 'andorra'), returnclass = "sf") %>%
  dplyr::select(name_en,
                iso_a2_eh,
                geometry) %>%
  rename(CNTR_CODE = iso_a2_eh,
         NAME_LATN = name_en) %>%
  mutate(NUTS_ID = CNTR_CODE)
