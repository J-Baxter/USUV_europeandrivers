####################  dependencies #################### 
library(ape)
library(tidyverse)

#################### import alignments, data and trees #################### 
# alignments
aln_files <- c('./2024Aug13/alignments/USUV_2024Aug13_alldata_aligned_formatted_noFLI_nflg.fasta',
               './2024Aug13/alignments/USUV_2024Aug13_alldata_aligned_formatted_noFLI_ns5_9100-9600.fasta')

aln_names <- gsub('.*noFLI_|\\.fasta', '' , aln_files)

aln <- lapply(aln_files,
              read.dna,
              format = 'fasta', 
              as.matrix = T)

#trees
tree_files <- c('./2024Aug13/iqtree_alldata/USUV_2024Aug13_alldata_aligned_formatted_noFLI_nflg.fasta.treefile',
                './2024Aug13/iqtree_alldata/USUV_2024Aug13_alldata_aligned_formatted_noFLI_ns5_9100-9600.fasta.treefile')

trees <- lapply(tree_files,
                read.tree)

#data
data <- read_csv('./data/metadata_noFLI_2024Aug13.csv')


#################### split data by alignment #################### 
data_per_alignment <- lapply(aln, 
                             function(x) data %>% 
                               filter(tipnames %in% rownames(x))) %>%
  setNames(aln_names) %>%
  bind_rows(, .id = 'alignment')

# sanity check
stopifnot(all(unlist(lapply(aln, nrow)) == unlist(lapply(data_per_alignment, nrow))))


#################### Subsampling 1 #################### 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across the entire 'global' dataset.

ExtractTerminalBranchLengths <- function(tree, aln_length, snp_threshold = 0){
  
  #create tidy tree object
  tidytree <- as_tibble(tree)
  
  sequence_groups <- tidytree %>%
    mutate(snps = floor(branch.length * aln_length)) %>% 
    filter(snps <= snp_threshold & grepl('\\|', label)) %>% 
    group_by(parent) %>%
    mutate(n = n()) %>% 
    filter(n > 1) %>%
    mutate(sequence_group = cur_group_id()) %>%
    ungroup() %>%
    select(c(label, sequence_group)) %>%
    rename(tipnames = label)
    
    
  return(sequence_groups)
}

groups <- mapply(ExtractTerminalBranchLengths,
                 trees,
                 lapply(aln, ncol),
                 snp_threshold = 0,
                 SIMPLIFY = F) %>%
  setNames(aln_names) %>%
  bind_rows(., .id = 'alignment')

data_per_alignment_wgroups <- data_per_alignment %>%
  left_join(groups)

subsample_1 <- data_per_alignment_wgroups %>%
  group_by(alignment) %>%
  
  # infer best location code and create grouping factor
  mutate(best_location_code = coalesce(collection_subdiv1code,
                                       collection_countrycode)) %>%
  mutate(sequence_group = as.factor(sequence_group)) %>%
  
  # group by sequence_identity > best_location > host_class
  
  group_by(alignment,
           sequence_group,
           best_location_code,
           host_class) %>%
  
  slice(which(date_ym == max(date_ym) | date_ym == min(date_ym))) %>%
  distinct() %>%
  ungroup() %>%
  group_split(alignment) %>%
  setNames(aln_names)

#nflg -> 126 sequences
#ns5_9100-9600 -> 302 sequences
# NB some identical sequences preserved (min and max within same location code)


#################### Subsampling 2 #################### 
# Details subsampling with respect to region. Here, the ML tree is traversed such that each time 
# region changes, the two most distant (temporally) tips are selected until the region next changes.
# This is applied across the entire 'global' dataset.

test <- as_tibble(t) %>%
  rename(tipnames = label) %>%
  left_join(data_per_alignment_wgroups %>% 
              filter(alignment == 'nflg') %>% 
              select(c(tipnames, collection_regionname))) %>%
  filter(grepl('\\|', tipnames)) %>%
  mutate(collection_regionchange = cumsum(collection_regionname != lag(collection_regionname, def = first(collection_regionname)))) %>% 
  group_by(Tag, LocationChange, location)
  


subsample_1 <- data_per_alignment_wgroups %>%
  #(left join tidytree with tr)
  
  group_by(alignment) %>%
  
  # infer best location code and create grouping factor
  mutate(best_location_code = coalesce(collection_subdiv1code,
                                       collection_countrycode)) %>%
  mutate(sequence_group = as.factor(sequence_group)) %>%
  
  # group by sequence_identity > best_location > host_class
  
  group_by(alignment,
           sequence_group,
           best_location_code,
           host_class) %>%
  
  slice(which(date_ym == max(date_ym) | date_ym == min(date_ym))) %>%
  distinct() %>%
  ungroup() %>%
  group_split(alignment) %>%
  setNames(aln_names)

