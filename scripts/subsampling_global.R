####################  dependencies #################### 
library(ape)
library(tidyverse)

#################### import alignments, data and trees #################### 
# alignments
aln_files <- c('./2024Aug13/alignments/USUV_2024Aug13_alldata_aligned_formatted_noFLI_nflg.fasta',
               './2024Aug13/alignments/USUV_2024Aug13_alldata_aligned_formatted_noFLI_ns5_9100-9600.fasta')

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
                               filter(tipnames %in% rownames(x)))

# sanity check
stopifnot(all(unlist(lapply(aln, nrow)) == unlist(lapply(data_per_alignment, nrow))))


#################### Subsampling 1 #################### 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across the entire 'global' dataset.


#################### Subsampling 2 #################### 
# Details subsampling with respect to region. Here, the ML tree is traversed such that each time 
# region changes, the two most distant (temporally) tips are selected until the region next changes.
# This is applied across the entire 'global' dataset.

