library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)

metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')
nflg_ml <- read.newick('./2024Oct20/alignments/nflg_rooted')
  
test <- tibble(node=1:Ntip(nflg_ml), tipnames = nflg_ml$tip.label) %>%
  mutate(tipnames = gsub("'", '', tipnames))

included_meta <- metadata %>% 
  select(nuts0_name,
         lineage,
         sublineage,
         tipnames) %>%
  unite(joint_lineage, lineage,sublineage) %>%
  mutate(joint_lineage =ifelse(joint_lineage == 'NA_NA', NA_character_, joint_lineage)) %>%
  left_join(test, by = join_by(tipnames)) %>%
  drop_na(node)

p1 <- p + geom_point2(aes(subset=node==16), color='darkgreen', size=5)
p2 <- rotate(p1, 16)
flip(p2, 17, 21)


nflg_ml %>%
  full_join(included_meta) %>%
  ggtree() +
  #geom_point2(aes(subset=node==364), color='darkgreen', size = 5) +
  geom_tippoint(aes(colour = joint_lineage)) +
  #geom_hilight(node=364, fill="steelblue", width = 1) +
  theme_tree()


read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta.treefile') %>%
  drop.tip(c('AY453411|NA|NA|AT|NA')) %>%
  write.tree('./2024Oct20/alignments/nflg_edited.treefile')
           
           