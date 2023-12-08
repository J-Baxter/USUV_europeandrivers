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
  filter(seq.length >= 300) 


plt1 <-data %>%
  group_by(year= floor_date(collection.date.formatted, "year"), iso.country, is.europe) %>%
  summarise(n = n()) %>%
  ggplot() + 
  geom_tile(aes(x = year, y = iso.country, fill=n)) + 
  scale_x_date('Collection Year') + 
  scale_y_discrete('Country') + 
  scale_fill_distiller(palette = 'OrRd', direction = 1) + 
  my_theme + 
  facet_grid(rows = vars(is.europe), scales = 'free_y', space = 'free_y',switch="both") + 
  theme(legend.position = 'bottom', strip.placement = 'outside')

ggsave(plot = plt1, filename = 'plt1.jpeg', device = jpeg,  width = 170, height = 120,  units = 'mm', dpi = 350)


plt2 <- data %>%
  ggplot()+
  #geom_rect(aes(xmin = 9000, xmax = Inf, ymin = 0, ymax = Inf), alpha = 0.007, fill = 'lightgrey')+
  geom_histogram(aes(x = seq.length, fill = segment), binwidth = 250) +
  geom_vline(xintercept  = 9000, linetype = 'dashed') + 
  geom_vline(xintercept  = 300, linetype = 'dashed') + 
  annotate('text', label = 'NFLG', x= 9500, y = 350, size = 3) +
  scale_x_continuous('Sequence Length', expand = c(0,0), breaks = seq(0,10000, by= 1000)) + 
  scale_y_continuous('Count', expand = c(0,0)) +
  scale_fill_brewer(palette = 'PuBu', na.value = 'black', 'Genomic Region') + 
  my_theme + 
  theme(legend.position = 'bottom')
  
ggsave(plot = plt2, filename = 'plt2.jpeg', device = jpeg,  width = 170, height = 120,  units = 'mm', dpi = 350)





  #group_by(month = floor_date(collection.date.formatted, "season"), iso.country) %>%
 # summarise(n = n()) %>%
  #ggplot() +
  #geom_line(aes(x = month, y = n, colour = iso.country))



#plt2 <- data %>%  
  #filter(seq.length >= 300) %>%#
  #filter(collection.date.formatted > '1990-01-01') %>%
  #group_by(host.order) %>%
  #filter(n() >= 5) %>%
  #ungroup() %>%
 # ggplot() +
  #geom_histogram(aes(x = collection.date.formatted, fill = host.order)) +
  #scale_fill_brewer(palette = 'BuGn', na.value = 'black') + 
  #facet_wrap(.~iso.country, scales = 'free_y') + 
  #my_theme + 
  #theme(legend.position = 'bottom')

#ggsave(plot = plt2, filename = 'plt2.jpeg', device = jpeg,  width = 190, height = 150,  units = 'mm', dpi = 350)


plt3 <- data %>%  
  filter(seq.length >= 300) %>%
  filter(host.class == 'Aves') %>%
  group_by(host.order) %>%
  summarise(n = n()) %>%
  mutate(host.freq = n / sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = host.order, y = host.freq, fill = '#fdbb84'), stat = 'identity') +
  scale_x_discrete('Host Order') + 
  scale_y_continuous(expand = c(0,0), 'Relative Frequency', limits =c(0,1.05), breaks = seq(0,1, by = 0.2))+
  my_theme 

ggsave(plot = plt3, filename = 'plt3.jpeg', device = jpeg,  width = 150, height = 150,  units = 'mm', dpi = 350)

plt4 <- data %>%  
  filter(host.order == 'Diptera') %>%
  group_by(host.genus, host.species) %>%
  #filter(n() > 5) %>%
  #group_by(host.species) %>%
  summarise(n = n()) %>%
  ungroup(host.genus) %>%
  mutate(host.freq = n / sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = host.species, y = host.freq, fill = '#fdbb84'), stat = 'identity') +
  scale_x_discrete('Mosquito Species') + 
  scale_y_continuous(expand = c(0,0), 'Relative Frequency', limits =c(0,1.05), breaks = seq(0,1, by = 0.2))+
  #scale_fill_brewer(palette = 'Paired') + 
  facet_grid(cols = vars(host.genus), scales = 'free_x', space = 'free_x', switch="both") +
  my_theme + 
  theme( strip.placement = 'outside') 



ggsave(plot = cowplot::plot_grid(plt3, plt4, align = 'hv', labels = c('Bird Order', 'Vector Species'), label_size = 10, label_fontfamily = "LM Sans 10", nrow = 2, hjust = -0.9), filename = 'plt4.jpeg', device = jpeg,  width = 190, height = 220,  units = 'mm', dpi = 350)