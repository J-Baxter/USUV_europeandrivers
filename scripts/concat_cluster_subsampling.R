####################################################################################################
####################################################################################################
## Script name: Subsample concatenated alignments (within-europe clusters)
##
## Purpose of script:
##
## Date created: 2024-10-31
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)
library(igraph)
library(phangorn)
library(phytools)

# User functions
# User functions
# Calculates pairwise hamming distance, then groups at a threshold difference
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


############################################## DATA ################################################
concat_alignment_files <- c('./2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_concat_A.fasta',
                          './2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_concat_B.fasta',
                          './2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_concat_C.fasta',
                          './2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI_concat_D.fasta')

concat_subset_alignments <- lapply(concat_alignment_files,
                                 read.dna,
                                 as.matrix = T,
                                 format = 'fasta')



# import metadata
concat_metadata <- read_csv('./data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv')

metadata_split <- lapply(concat_subset_alignments, function(x) metadata %>% filter(tipnames %in% rownames(x)))

############################################## MAIN ################################################

# exclude Tempest outliers

to_exclude <- c(
  # A - unclocklike
  'JF834629|USU-MO-m9/2010|culex_pipiens|IT|2010-08-19',
  'KY113104|BNI-527|NA|DE|2016',
  'MK598468|NA|turdidae_sp|AT|2018-08-16',
  'KY113101|BNI-504|NA|DE|2016',
  'KY113103|BNI-511|NA|DE|2016',
  'HM138709|USU-MO-1/2009|culex_pipiens|IT|2009-09-09',
  'MK598463|NA|turdus_merula|AT|2018-08-08',
  'MW803140|UV780/Serbia/2019|culex_pipiens|RS|2019',
  'MN395372|NA|turdus_merula|CZ|NA',
  
  # B - unlocklike
  'MK230893|Grivegnee/2017|turdus_merula|BE|2017-08',
  'NA|USUV-6736/France/2017|strix_nebulosa|FR|2017-08',
  'NA|USUV-5442/France/2018|tetrao_urogallus|FR|2018-08-21',
  'NA|USUV-10212/France/2016|turdus_merula|FR|2016-09-22',
  'MK060108|NA|homo_sapiens|IT|2017-08',
  
  
  
  # C unclocklike
  'MG461314|USUV_379m_Israel04_USUTU_|culex_pipiens|IL|2004-10-14',
  'NA|USUV-8784/France/2016|turdus_merula|FR|2016-09-23',
  'NA|USUV-9784/France/2017|turdidae_sp|FR|2017-09-30',
  'NA|USUV-7391/France/2017|turdus_merula|FR|2017-09-12',
  'NA|USUV-7641/France/2017|turdus_merula|FR|2017-09-18',
  'NA|USUV-8760/France/2016|turdus_merula|FR|2016-09-21')

concat_subset_alignments %<>% 
  lapply(., function(x) x[!rownames(x) %in% to_exclude,])


# Subsampling 1 
# Details stratified subsampling aiming to 'thin' the tree by removing sequences that are identical
# in composition and near-identical in traits.
# This is applied across each cluster alignment

groups <- lapply(concat_subset_alignments,
                 GroupSequences,
                 snp_threshold = 1) %>%
  
  setNames(., LETTERS[1:4]) %>%
  bind_rows(., .id = 'CLUSTER_GROUP')



subsample_1 <- metadata_split %>%
  
  # set names
  setNames(., LETTERS[1:4]) %>%
  
  bind_rows(., .id = 'CLUSTER_GROUP') %>%
  
  # join identity groups
  left_join(groups) %>%
  
  # establish groupings
  group_by(CLUSTER_GROUP,
           nuts0_id,
           date_y,
           date_ym,
           nuts1_id) %>%
  
  #slice(which(date_y == max(date_y, na.rm = T) | date_y == min(date_y, na.rm = T)), 
       # .preserve = T) %>%
 # group_by(date_y, .add = TRUE) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  group_split(CLUSTER_GROUP)


aln_subsampled <- mapply(function(aln, meta) aln[rownames(aln) %in% meta$tipnames,],
                         concat_subset_alignments,
                         subsample_1)




beast_data <- subsample_1 %>% 
  
  mapply(function(aln, meta) meta %>% filter(tipnames %in% rownames(aln)),
         concat_subset_alignments,
         .,
         SIMPLIFY = F) %>%
  
  # set names
  setNames(., LETTERS[1:4]) %>%
  do.call(rbind.data.frame,. ) %>%

  
  # spread geocode coords to lat lon
  mutate(geocode_coords = gsub('c|\\(|\\)|,', '', geocode_coords)) %>%
  separate_wider_delim(geocode_coords, ' ', names =  c('lon', 'lat')) %>%
  mutate(across(c(lat, lon), .fns = ~ as.numeric(.x))) %>%
  
  
  # Select columns for BEAST
  dplyr::select(
    tipnames,
    lat,
    lon,
    nuts0_id,
    nuts1_id,
    CLUSTER_GROUP) %>%
  group_split(CLUSTER_GROUP, .keep = FALSE)



############################################## WRITE ###############################################



############################################## WRITE ###############################################
filenames <- gsub('.fasta', '_subsampled.fasta', concat_alignment_files)

mapply(write.dna,
       aln_subsampled,
       filenames,
       format = 'fasta')

beast_filename <- gsub('_alldata_aligned_formatted', '', concat_alignment_files) %>%
  gsub('.fasta$', '_traits.txt', .)

mapply(write_delim,
       beast_data,
       beast_filename,
       delim = '\t',
       quote= 'needed')
# A
#MN419903|NA|turdus_merula|CZ|2018-07-16
#MK060111|NA|culex_sp|IT|2018-09
#HM138709|USU-MO-1/2009|culex_pipiens|IT|2009-09-09  
#HM138716|USU-BO-8/2009|culex_pipiens|IT|2009-09-02  
# C
#MG461314|USUV_379m_Israel04_USUTU_|culex_pipiens|IL|2004-10-14    
#ON758918|IREC01|turdus_merula|ES|2018                             
#ON758919|IREC02|turdus_merula|ES|2018                            
#ON099437|Usutu/Camargue2018|culex_pipiens|FR|2018                 
#ON099439|Usutu/Camargue2020/poolB5|culex_pipiens|FR|2020          
#ON099438|Usutu/Camargue2020/poolH4|culex_pipiens|FR|2020          
#MH423836|DoMaPS455/2016|culex_pipiens|DE|2016-07                
#NA|USUV-6056/France/2018|turdus_merula|FR|2018-09-05  

# B
#MN419907|NA|turdus_merula|CZ|2018-09-18                                                 
#MN419909|NA|turdus_merula|CZ|2018-09-11                                                   
#MN419910|NA|turdus_merula|CZ|2018-09-04                                                  
#MN419911|NA|turdus_merula|CZ|2018-09-05                                                  
#MN419912|NA|turdus_merula|CZ|2018-09-17                                                   
#MN419908|NA|turdus_merula|CZ|2018-09-18                                                   
#KU664615|27-Nov|strix_nebulosa|DE|2011                                                    
#KU664612|24-Nov|strix_nebulosa|DE|2011                                                    
#KY114797|NA|strix_nebulosa|DE|2016                                                       
#KU664610|21-Nov|strix_nebulosa|DE|2011                                                    
#KU664611|22-Nov|strix_nebulosa|DE|2011                                                   
#KU664614|26-Nov|strix_nebulosa|DE|2011                                                    
#KU664613|25-Nov|strix_nebulosa|DE|2011 
# D
#ON125263|Usutu/Camargue/Culex2015|culex_pipiens|FR|2015  
#KU664609|Usutu_Berlin|strix_nebulosa|DE|2015-08-27       