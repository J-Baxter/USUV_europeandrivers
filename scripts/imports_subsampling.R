#nflg downsampling

#dependencies
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(ggtree)


# import ML tree
nflg_ml <- read.newick('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta.contree')

# import metadata
metadata <- read_csv('./data/USUV_metadata_all_2024Oct20.csv')


############# Import 'Master' Alignment ############# 
nflg_alignment <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta',
                as.matrix = T,
                format = 'fasta')

# exclude Tempest outliers
to_exclude <- c('MK230893|Grivegnee/2017|turdus_merula|BE|2017-08',
                'OQ414983|TV193/2019|culex_pipiens|RO|2019-07',
                'KT445930|13-662|culex_modestus|CZ|2013-08-07')

nflg_alignment %<>% .[!rownames(.) %in% to_exclude,]



#################### Subsampling 1 #################### 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across the entire 'global' dataset.

# Calculates pairwise hamming distance, then groups at a threshold difference
GroupSequences <- function(aln, snp_threshold = 0){
  
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
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}


groups <- GroupSequences(nflg_alignment, 5) 

groups %>%
  summarise(n = n(), .by = sequence_group) %>%
  pull(n) %>%
  hist()


subsample_1 <- metadata %>%
  
  # include only sequences in alignment
  filter(tipnames %in% rownames(nflg_alignment)) %>%
  
  # join identity groups
  left_join(groups) %>%
  
  # establish groupings
  group_by(nuts1_id,
           date_y,
           sequence_group) %>%
  
  slice_sample(n = 1)


aln_subsampled <- nflg_alignment[rownames(nflg_alignment) %in% subsample_1$tipnames,]

write.dna(aln_subsampled,
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_nflg_subsample1.fasta', 
          format = 'fasta')



#################### Subsampling 2 #################### 

