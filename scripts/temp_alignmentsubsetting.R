FindIdenticalSeqs<- function(alignment, snp_threshold = 1){
  
  # calculate pairwise distances
  distances <- ape::dist.dna(alignment, model = "N")
  
  # Identify identical sequences
  identical_pairs <- which(as.matrix(distances) < 3, arr.ind = TRUE) %>%
    dplyr::as_tibble(rownames = 'tipnames') %>% 
    filter(row !=col) 
  
  # infer adjcency network 
  sequence_net <- igraph::graph_from_data_frame(identical_pairs %>% 
                                                  select(!tipnames), 
                                                directed = F)
  
  # returns a dataframe of tipnames and group number (int)
  groupings <- merge(identical_pairs,
                     stack(igraph::clusters(sequence_net)$membership),
                     by.x = "row", 
                     by.y = "ind", all.x = TRUE) %>%
    dplyr::select(c(tipnames, values)) %>%
    dplyr::distinct() %>%
    dplyr::rename(group = values)
  
  return(groupings)
}

alignments = list.files()
nflg <- read.dna('./2024Apr21/alignments/USUV_nflg_2024May7.fasta',
                 format = 'fasta',
                 as.matrix= T)

identical_nflg <- FindIdenticalSeqs(nflg)

nflg_data_subset <-  data_formatted %>% left_join(., identical_nflg) %>% slice_sample(n=1, by = c(group, collection_countrycode))

nflg_subsampled_alignment <- nflg [rownames(nflg ) %in% nflg_data_subset$tipnames,]


write.FASTA(nflg_subsampled_alignment,
            './2024Apr21/alignments/USUV_nflg_2024May7_subsampled.fasta'
            )



e <- read.dna('./2024Apr21/alignments/USUV_E_2024May7.fasta',
                 format = 'fasta',
                 as.matrix= T)

identical_e <- FindIdenticalSeqs(e)

e_data_subset <-  data_formatted %>% left_join(., identical_nflg) %>% slice_sample(n=1, by = c(group, collection_countrycode))

e_subsampled_alignment <- e[rownames(e) %in% e_data_subset$tipnames,]


write.FASTA(e_subsampled_alignment,
            './2024Apr21/alignments/USUV_e_2024May7_subsampled.fasta'
)