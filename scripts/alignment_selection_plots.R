# read_alignments
library(ape)
aln_files <- list.files('./2024Aug13/alignments/',
                        full.names = T,
                        pattern = 'fasta')[1:5]

aln <- lapply(aln_files, read.dna, format = 'fasta', as.matrix = TRUE)

data <- read_csv('./data/metadata_noFLI_2024Aug13.csv')

data_per_alignment <- lapply(aln, function(x) data %>% filter(tipnames %in% rownames(x))) %>%
  setNames(gsub('.*USUV_2024Aug13_alldata_aligned_formatted_noFLI_|\\.fasta',  '', aln_files)) %>%
  bind_rows(, .id = 'aln')

alignmentinfo <- data_per_alignment %>%
  group_by(aln, ) %>%
  summarise(n = n(), n_countries = n_distinct(collection_countryname), earliest_year = min(date_y), mostrecent_year = max(date_y), list_countries = paste(unique(collection_countryname), collapse = ','))


data_per_alignment %>%
  group_by(aln, year= floor_date(date_ymd, "year"), collection_countryname) %>%
  summarise(n = n()) %>%
  ggplot() + 
  geom_tile(aes(x = year, y = collection_countryname, fill=n)) + 
  scale_x_date('Collection Year') + 
  scale_y_discrete('Country') + 
  scale_fill_viridis_c(option = 'B') +
  #scale_fill_distiller(palette = 'OrRd', direction = 1) + 
  #my_theme + 
  facet_grid(rows = vars(aln), scales = 'free', space = 'free',switch="both", drop = T) + 
  theme(legend.position = 'bottom', strip.placement = 'outside')