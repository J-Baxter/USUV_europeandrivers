#### Classify partial sequences according to old lineages 
# uses molecular data published by Zecchin et al. and Engel et al. 


# NFLGY ONLY
lineage_signatures <- read_csv('./data/zechchin_lineagesignatures.csv') %>%
  mutate(across(-1, .fns = ~ gsub(',', '|', .x)))


############# Import 'Master' Alignment and metadata ############# 
aa_aln <- read.dna('./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta',
                as.matrix = T,
                format = 'fasta') %>%
  trans() %>%
  as.character.AAbin()


metadata_noFLI <- read_csv('./data/USUV_metadata_noFLI_2024Oct20.csv')

GetLineage <- function(seq, lineage_signatures){

  lineage <- apply(lineage_signatures[,-1],
                   2,
                   function(x) mapply(grepl,x, seq[lineage_signatures$codon])) %>%
    apply(., 2, all) %>%
    which() %>%
    names()
  
  if(is.null(lineage)){
    lineage <- NA
  }
  
  return(lineage)
}


#key_sites <- aa_aln[, lineage_signatures$codon]

test_lineages <- apply(aa_aln, 1, GetLineage, lineage_signatures) %>%  
  as.matrix() %>%
  as_tibble(rownames = 'tipnames') %>%
  mutate(V1 = as.character(V1) %>%
           gsub('character\\(0\\)', NA_character_,.)) %>%
  rename(lineages = V1)


test <- metadata %>% dplyr::select(c(tipnames,
                                     lineage,
                                     sublineage,
                                     generegion_nflg)) %>%
  filter(generegion_nflg == 1) %>%
  left_join(test_lineages)


#Stopped because protein sequences are not distinguishable between lineages in short NS5 segments

