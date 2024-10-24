############# Dependencies ############# 
library(tidyverse)
library(ape)


############# Import 'Master' Alignment ############# 
aln <- read.dna('./2024Aug13/alignments/USUV_2024Aug28_alldata_aligned_formatted_noFLI.fasta',
                as.matrix = T,
                format = 'fasta')



############# Extract Sub Alignment ############# 
sub_alignment_coords <- tibble(start = c(1003, 8968, 9100, 10042),
                               region = c('env', 'ns5', 'ns5', 'ns5'),
                               end = c(1491, 9264, 9600, 10314))

sub_alignments <- apply(sub_alignment_coords, 
                        1, 
                        function(x) aln[,x['start']:x['end']]) %>%
  lapply()


############# Write Sub Alignment to File ############# 

sub_alignment_names <- paste0(
  './2024Aug13/alignments/USUV_2024Aug28_alldata_aligned_formatted_noFLI_',
  sub_alignment_coords$region, '_',
  sub_alignment_coords$start, '-',
  sub_alignment_coords$end,
  '.fasta')

mapply(write.dna,
       sub_alignments,
       sub_alignment_names,
       format = 'fasta')



###

#############  Sub-aligments ############# 
write.dna(alignment[rownames(alignment) %in% (metadata_noFLI %>% 
                                                filter(generegion_nflg == 1) %>% 
                                                pull(tipnames)),],
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NFLG.fasta', 
          format = 'fasta')


# Alignment has 953 sequences with 665 columns, 385 distinct patterns, 152 parsimony-informative, 80 singleton sites, 433 constant sites
write.dna(alignment[rownames(alignment) %in% (metadata_noFLI %>% 
                                                filter(generegion_NS5_9000_9600 == 1) %>% 
                                                pull(tipnames)),],
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NS5_9000-9600.fasta', 
          format = 'fasta')

# 
#write.dna(alignment[rownames(alignment) %in% (metadata_noFLI %>% 
#  filter(generegion_NS5_9100_9600 == 1) %>%
# pull(tipnames)),],
#'./2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NS5_9100-9600.fasta', 
# format = 'fasta')


write.dna(alignment[rownames(alignment) %in% (metadata_noFLI %>% 
                                                filter(generegion_NS5_10042_10312 == 1) %>% 
                                                pull(tipnames)),],
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_NS5_10042-10312.fasta', 
          format = 'fasta')


write.dna(alignment[rownames(alignment) %in% (metadata_noFLI %>%
                                                filter(generegion_env_1003_1491 == 1) %>%
                                                pull(tipnames)),],
          './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_env_1003-1491.fasta', 
          format = 'fasta')
