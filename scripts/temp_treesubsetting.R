library(ggtree)
library(ggtreeExtra)
library(treeio)
library(TreeTools)
my_theme <- theme_classic()+ #base_family = "LM Sans 10"
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

offspring.tbl_tree_item <- getFromNamespace("offspring", "tidytree")
assign("offspring.tbl_tree_item", offspring.tbl_tree_item, envir = .GlobalEnv)
known_lineages <- read_csv('./data/known_lineages.csv')

nflg_tree <- read.beast('./2024Apr21/USUV_nflg_2024May7_subsampled_mcc') %>% left_join(., data_formatted, by = join_by(label == tipnames))%>%
  left_join(., 
            known_lineages %>% filter(!is.na(accession)), 
            by = join_by( sequence_accession == accession), relationship = 'one-to-one') %>%
  left_join(., known_lineages %>% filter(!is.na(isolate)), 
            by = join_by( sequence_isolate == isolate)) %>%
  mutate(lineage = coalesce(lineage, lineage.y))

nflg_data <- data.frame(
  node = c(166, 197, 215, 141, 125, 118), 
  name = c("EU3", "EU2", 'EU1', 'AF3.3', 'AF3.2', 'AF3.1'))



ggtree(nflg_tree, mrsd="2023-09-23", right = TRUE) + 
  theme_tree2() +
  geom_tippoint(aes( color=collection_countryname))  +
  geom_tiplab(aes(label = lineage)) +
  geom_nodelab(aes(label = node)) +
  geom_cladelab(data = nflg_data, 
                  mapping = aes(node = node, label = name), size=18,
                  align = T) + 
  theme(axis.text =  element_text(size=18)) 


eu3 <- as_tibble(nflg_tree) %>% 
  offspring(166) %>%
  pull(label) %>% 
  .[!is.na(.)]

eu3_groups <- metadata %>% 
  left_join(., identical_nflg) %>%
  filter(tipnames %in% eu3) %>%
  pull(group) %>%
  unique()


eu3_data_sample <-  metadata %>% 
  left_join(., identical_nflg) %>%
  filter(group %in% eu3_groups) %>%
  mutate(collection_bestlocation = coalesce(collection_subdiv2code, collection_subdiv1code, collection_countrycode)) %>%
  #filter(!is.na(collection_subdiv1name)) %>%
  slice_sample(n=1, by = c(group, collection_bestlocation))

EU3_nflg_subsampled_alignment <- nflg [rownames(nflg ) %in% eu3_data_sample$tipnames,]


write.FASTA(EU3_nflg_subsampled_alignment,
            './2024Apr21/alignments/USUV_EU3_nflg_2024May7_subsampled.fasta'
)
  
t <- eu3_data_sample %>%
  select(c(tipnames, starts_with('collection'))) %>%
  mutate(collection_subdiv2long = case_when(collection_subdiv2long == 'list' ~ NA, .default = collection_subdiv2long)) %>%
  mutate(collection_bestlat = coalesce(collection_subdiv2lat, collection_subdiv1lat, collection_countrylat)) %>%
  mutate(collection_bestlong = coalesce(collection_subdiv2long, collection_subdiv1long, collection_countrylong)) %>%
  select(c(tipnames, starts_with('collection_best'))) %>%
  mutate(across(c(ends_with('lat'), ends_with('long')), .fns = ~as.double(.x) )) %>%

  write_delim(delim = '\t', quote= 'needed',
             'USUV_EU3_nflg_2024May7_subsampled.txt')    


e_tree <- read.beast('./2024Apr21/USUV_E_2024May7_subsampled_mcc') %>% left_join(., data_formatted, by = join_by(label == tipnames)) %>%
  left_join(., 
            known_lineages %>% filter(!is.na(accession)), 
            by = join_by( sequence_accession == accession), relationship = 'one-to-one') %>%
  left_join(., known_lineages %>% filter(!is.na(isolate)), 
            by = join_by( sequence_isolate == isolate)) %>%
  mutate(lineage = coalesce(lineage, lineage.y))
  
e_data <- data.frame(
  node = c(166, 197, 215, 141, 125, 118), 
  name = c("EU3", "EU2", 'EU1', 'AF3.3', 'AF3.2', 'AF3.1'))


ggtree(e_tree, mrsd="2023-09-22", right = TRUE) + 
  theme_tree2() + 
  geom_tippoint(aes( color=collection_countryname))  +
  geom_tiplab(aes(label = lineage)) +
  geom_cladelab(data = e_data, 
                mapping = aes(node = node, label = name), 
                size=18,
                align = T) + 
  theme(axis.text =  element_text(size=18)) 

ggtree(Subtree(e_tree@phylo, 117), mrsd="2023-09-22", right = TRUE) + 
  theme_tree2()


geom_cladelabel()


test <- data_formatted %>%
  #left_join(., known_lineages, 
         #   by = join_by( sequence_isolate == isolate)) %>%
  left_join(., known_lineages, 
            by = join_by( sequence_accession == accession)) 
  


t <- nflg_tree@extraInfo
