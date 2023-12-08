# USUV public data analysis

library(tidyverse)
library(scales)
library(extrafont)


#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 
loadfonts()

my_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )

data <- read_csv('./data/genbank_usuv_20231206.csv') %>%
  select(-1) %>% 
  filter(seq.length >= 300) %>%
  filter(collection.date.formatted > '1980-01-01') %>%
  group_by(host.order) %>%
  filter(n() >= 5) %>%
  ungroup() 
  

data %>%
  group_by(month = floor_date(collection.date.formatted, "season"), iso.country) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(x = month, y = n, colour = iso.country))



plt2 <- data %>%  
  filter(seq.length >= 300) %>%
  filter(collection.date.formatted > '1990-01-01') %>%
  group_by(host.order) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  ggplot() +
  geom_histogram(aes(x = collection.date.formatted, fill = host.order)) +
  scale_fill_brewer(palette = 'BuGn', na.value = 'black') + 
  facet_wrap(.~iso.country, scales = 'free_y') + 
  my_theme + 
  theme(legend.position = 'bottom')

ggsave(plot = plt2, filename = 'plt2.jpeg', device = jpeg,  width = 190, height = 150,  units = 'mm', dpi = 350)

plt3 <- data %>%  
  filter(seq.length >= 300) %>%
  filter(collection.date.formatted > '1990-01-01') %>%
  filter(host.class == 'Aves') %>%
  group_by(host.order) %>%
  filter(n() > 5) %>%
  group_by(host.genus) %>%
  summarise(n = n()) %>%
  mutate(host.freq = n / sum(n)) %>%
  ungroup() %>%
  left_join(data %>% select(c(host.genus, host.order)) %>% distinct()) %>%
  ggplot() +
  geom_bar(aes(x = host.order, y = host.freq, fill = host.genus), stat = 'identity') +
  #scale_fill_brewer(palette = 'Paired', na.value = 'black') + 
  #facet_wrap(.~iso.country, scales = 'free_y') + 
  my_theme + 
  theme(legend.position = 'right')

ggsave(plot = plt3, filename = 'plt3.jpeg', device = jpeg,  width = 150, height = 150,  units = 'mm', dpi = 350)

plt4 <- data %>%  
  filter(seq.length >= 300) %>%
  filter(collection.date.formatted > '1990-01-01') %>%
  filter(host.order == 'Diptera') %>%
  group_by(host.genus) %>%
  filter(n() > 5) %>%
  group_by(host.species) %>%
  summarise(n = n()) %>%
  mutate(host.freq = n / sum(n)) %>%
  ungroup() %>%
  inner_join(data %>% 
             filter(host.order == 'Diptera') %>%
               group_by(host.genus) %>%
               filter(n() > 5) %>% 
               ungroup() %>% 
               select(c(host.species, host.genus)) %>%
               distinct() %>% 
               filter(if_any(everything(), ~ !is.na(.))))  %>%
  ggplot() +
  geom_bar(aes(x = host.genus, y = host.freq, fill = host.species), stat = 'identity') +
  scale_fill_brewer(palette = 'BuGn', na.value = 'black') + 
  #facet_wrap(.~iso.country, scales = 'free_y') + 
  my_theme + 
  theme(legend.position = 'right')

ggsave(plot = plt4, filename = 'plt4.jpeg', device = jpeg,  width = 150, height = 150,  units = 'mm', dpi = 350)