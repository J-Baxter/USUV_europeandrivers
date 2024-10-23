library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)

metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')

nflg_ml <- read.newick('./2024Oct20/alignments/nflg_rooted')


test <- tibble(node=1:Ntip(nflg_ml), tipnames = nflg_ml$tip.label) %>%
  mutate(tipnames = gsub("'", '', tipnames))

included_meta <- metadata %>% 
  dplyr::select(nuts0_name,
         lineage,
         sublineage,
         tipnames) %>%
  unite(joint_lineage, lineage,sublineage, remove = F) %>%
  mutate(joint_lineage =ifelse(joint_lineage == 'NA_NA', NA_character_, joint_lineage)) %>%
  left_join(test, by = join_by(tipnames)) %>%
  drop_na(node)


 nflg_ml %>%
  full_join(included_meta) %>%
  ggtree() +
  #geom_point2(aes(subset=node==364), color='darkgreen', size = 5) +
  geom_tippoint(aes(colour = lineage)) +
  #geom_hilight(node=364, fill="steelblue", width = 1) +
  theme_tree() + 
   geom_nodelab(size = 2)+
  geom_tiplab(size = 2) 


read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta.treefile') %>%
  drop.tip(c('AY453411|NA|NA|AT|NA')) %>%
  write.tree('./2024Oct20/alignments/nflg_edited.treefile')


ns5_ml <- read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NS5_9000-9600.fasta.contree')


ns5_ml %>%
  full_join(included_meta) %>%
  ggtree() +
  #geom_point2(aes(subset=node==364), color='darkgreen', size = 5) +
  geom_tippoint(aes(colour = lineage)) +
  #geom_hilight(node=364, fill="steelblue", width = 1) +
  theme_tree() + 
  geom_tiplab(size = 1) 

####### Subample 1 tree
subsample_ML <- read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta.treefile')

test <- tibble(node=1:Ntip(subsample_ML), tipnames = subsample_ML$tip.label) %>%
  mutate(tipnames = gsub("'", '', tipnames))

included_meta <- metadata %>% 
  dplyr::select(nuts0_id,
                is_europe,
                lineage,
                sublineage,
                tipnames) %>%
  unite(joint_lineage, lineage,sublineage, remove = F) %>%
  mutate(joint_lineage =ifelse(joint_lineage == 'NA_NA', NA_character_, joint_lineage)) %>%
  left_join(test, by = join_by(tipnames)) %>%
  drop_na(node)

subsample_ML %>%
  midpoint() %>%
  full_join(included_meta) %>%
  ggtree() +
  #geom_point2(aes(subset=node==364), color='darkgreen', size = 5) +
  geom_tippoint(aes(colour = as.factor(is_europe))) + 
  geom_cladelabel(node=581, 'AF 3 ', align = T) + #africa 3
  geom_cladelabel(node=520, 'EU 3 ', align = T) + # europe 3
  geom_cladelabel(node=574, 'EU 4 ', align = T) + #europe 4
  geom_cladelabel(node=516, 'EU 1 ', align = T) + #europe 1
  geom_cladelabel(node=652, 'AF 2 ', align = T) + #africa 2
  geom_cladelabel(node=501, 'EU 2 ', align = T) + #europe 2
  theme_tree() +
  #geom_text(aes(label=node)) + 
  theme(legend.position = 'bottom') +
  scale_colour_brewer(NULL, palette = 'Dark2', labels = c('0' = 'Outside of Europe',
                                                          '1' = 'Within Europe'))

           