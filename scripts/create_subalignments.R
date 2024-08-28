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