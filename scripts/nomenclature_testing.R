####################################################################################################
####################################################################################################
## Script name: Tree traversal + Nomenclature
##
## Purpose of script: Propose a suitabl nomenclautre based on Hill et al. 2023
##
## Date created: 2024-12-04
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions
GroupSequences <- function(aln, snp_threshold = 0){
  require(igraph)
  
  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln) 
  }else{
    aln_formatted <- aln
  }
  
  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  if( any(hd_raw[lower.tri(hd_raw, diag = FALSE)] <= snp_threshold)){
    groups <- which(hd_raw <= snp_threshold, 
                    arr.ind = TRUE) %>%
      dplyr::as_tibble(rownames = 'tipnames') %>% 
      filter(row !=col) %>%
      dplyr::select(-tipnames) %>%
      
      # Infer network from HDs
      igraph::graph_from_data_frame(.,
                                    directed = F) %>%
      
      components() %>%
      getElement('membership') %>%
      stack() %>%
      as_tibble() %>%
      mutate(ind = as.numeric(as.character(ind))) %>%
      mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
      dplyr::select(c(tipnames, values)) %>%
      dplyr::distinct() %>%
      dplyr::rename(sequence_group = values)
  }else{
    print('all sequences above threshold.')
    groups <- tibble(tipnames = rownames(aln)) %>%
      rowid_to_column(., var = 'sequence_group')
  }
  
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}


nflg_aln <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta',
                     format = 'fasta', 
                     as.matrix = T)



############################################## DATA ################################################
ml_tree <- read.tree('nextstrain/results/tree_raw.nwk')

#nflg_mcc <- read.beast('./2024Oct20/test_beast/USU')

metadata_in_tree <- read_csv('./data/USUV_metadata_all_2024Dec02.csv') %>%
  filter(tipnames %in%ml_tree$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(ml_tree)) #sanity check

############################################## MAIN ################################################
most_recent_date <- metadata_in_tree %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym - NB this no longer holds for hungary
  max(na.rm = TRUE)

tree_tbl <- as_tibble(ml_tree) %>%
  mutate(branchlength_subs = branch.length*10315)


# Pango-like lineages
#1. Minimum support of 70%
#2. Containing at least 5 sequences
#3. genomes within the lineage must share at least one nucleotide change relative to ancestral
ml_tree@phylo <- midpoint.root(ml_tree@phylo)
#http://127.0.0.1:19839/graphics/plot_zoom_png?width=1400&height=1190
wellsupported_nodes <- ml_tree %>% 
  as_tibble() %>%
  # Only include nodes with bootstrap support greater than 70%
  filter(UFboot >= 70) %>%
  pull(node)
  

#1. Minimum support of 70%
# n = 263
subtrees <- lapply(wellsupported_nodes, function(x) as_tibble(ml_tree) %>% offspring(., x, self_include = T)) %>%
  
#2. Containing at least 5 sequences
  # Force NA means we only count terminal nodes
# n = 119
  .[sapply(. , function(x) x %>% filter(is.na(UFboot)) %>% nrow(.) >= 5)] %>%
  
#3. Genomes within the lineage must share at least one nucleotide change relative to ancestral
  # (filter parent branch must be at least one substitution different)
# n = 118
  .[sapply(., function(x) x %>% 
           filter(!parent %in% node) %>% 
           pull(branch.length) %>%
           multiply_by(10315) %>%
             round() %>%
             is_weakly_greater_than(10))] %>%
  bind_rows(., .id = 'lineage')


lineage_nodes <- subtrees %>%
  group_by(lineage) %>%
  filter(!parent %in% node) %>% 
  ungroup() %>%
  select(node) %>%
  mutate(is_start_lineage = '1')


  

  
autolin_labels <- read_tsv('./2024Dec02/nomenclature/labels.tsv') %>%
  rename(label = sample) %>%
  separate_wider_delim(lineage, delim = '.', names = c('level_0', 'level_1', 'level_2', 'level_3'),
                       too_few = 'align_start') %>%
  mutate(level_3 = case_when(!is.na(level_3) ~ paste(level_0, level_1, level_2, level_3, sep = '.')),
         level_2 = case_when(!is.na(level_2) ~ paste(level_0, level_1, level_2, sep = '.')),
         level_1 = case_when(!is.na(level_1) ~ paste(level_0, level_1, sep = '.')))


midpoint.root(ml_tree) %>% 

  #left_join(metadata_in_tree %>% 
             # dplyr::select(is_europe, sequence_accession, nuts0_id, tipnames, lineage, host_class) %>%
            #  rename(label = tipnames),
          #  by = 'label') %>%
  
 # left_join(tree_tbl) %>% 
 # mutate(wellsupported = ifelse( UFboot > 0.7, '1', '0')) %>%
  
  left_join(autolin_labels) %>%
  
 # left_join(lineage_nodes) %>%
  ggtree() + 
  theme_tree2() + 
  
  geom_cladelabel()

  geom_fruit(geom = geom_tile,
             mapping = aes(fill = level_0),
             width = 0.002,
             colour = NA,
             pwidth = 1.2,
             offset = 0.03) + 
  geom_fruit(geom = geom_text,
             mapping = aes(label = level_0),
             offset = 0.03) + 
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = level_1),
             width = 0.002,
             colour = "white",
             pwidth = 1.2,
             offset = 0.06)
  scale_fill_brewer('Location', palette = 'Dark2', direction = -1, 
                    labels = c('0' = 'Outside of Europe',
                               '1' = 'Within Europe'),
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 2)) 
  #geom_text(aes(label=node), hjust=-.3)
  #scale_colour_distiller(palette = 'OrRd', direction = 1)
  
  
test <- treeio::read.nextstrain.json('./2024Dec02/nomenclature/annotated.json')

test_mrca_level0 <- test@data %>%
  filter(`GRI Lineage Level 0` != 'not assigned') %>%
  group_by(`GRI Lineage Level 0`) %>%
  slice_min(num_date) %>% 
  select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 0`)


test_mrca_level1 <- test@data %>%
  filter(`GRI Lineage Level 1` != 'not assigned') %>%
  group_by(`GRI Lineage Level 1`) %>%
  slice_min(num_date) %>% 
  select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 1`) %>%
  mutate(label = case_when(!grepl('\\.', label) ~ NA_character_, .default = label)) %>%
  drop_na()



test_mrca_level2 <- test@data %>%
  filter(`GRI Lineage Level 2` != 'not assigned') %>%
  group_by(`GRI Lineage Level 2`) %>%
  slice_min(num_date) %>% 
  select(node, country, starts_with('GRI')) %>%
  rename(label = `GRI Lineage Level 2`) %>%
  mutate(label = case_when(!grepl('\\.[[:digit:]]\\.', label) ~ NA_character_, .default = label)) %>%
  drop_na()


 test %>% 
  ggtree() + 
  geom_cladelab(data = test_mrca_level0 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 15,
                offset = 0.005,
                offset.text = -0.0004) + 
   geom_hilight(data = test_mrca_level0 , 
                mapping = aes(node = node, label = label, fill = label))+ 
  geom_cladelab(data = test_mrca_level1 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 15,
                offset = 0.003,
                offset.text = -0.0005) + 
  geom_cladelab(data = test_mrca_level2 , 
                mapping = aes(node = node, label = label, colour = label), 
                align = TRUE, 
                fontsize = 4, 
                textcolour = 'black',
                barsize = 15,
                offset = 0.001,
                offset.text = -0.0005)
  
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################