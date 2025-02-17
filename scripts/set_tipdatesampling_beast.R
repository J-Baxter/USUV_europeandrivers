# write nexus for beast2 analysis:
# requires charset
# requires calibrate for tipdate sampling
# requires additional txt file output for traits

library(ape)

aln <- read.dna('./2025Feb10/alignments/USUV_2025Feb10_alldata_aligned_formatted_noFLI_NFLG_subsample_p10.fasta',
         format = 'fasta',
         as.matrix = T)


WriteBEAST2Nexus <- function(aln,
                             filename,
                             metadata, 
                             SRD06 = TRUE){
  
  stopifnot(is.matrix(aln))
  
  # Extract key matrix info
  names <- rownames(aln)
  n_sequences <- nrow(aln)
  n_positions <- ncol(aln)
  
  metadata_in_aln <- metadata %>%
    filter(tipnames %in% names)
  
  # Sanity check
  stopifnot(nrow(metadata_in_aln) == n_sequences)

  
  
  # Specify header block
  header_block <- list('#NEXUS', 
                       '')
  
  # Format alignment block
  sequences_as_characters <- split(as.character(aln), row(as.character(aln))) %>%
    lapply(., paste, collapse = '') %>%
    flatten_chr() %>%
    str_to_upper()
  
  labelled_alignment <- list(names,  sequences_as_characters) %>%
    do.call(Map, c(c, .)) %>%
    flatten_chr()

  alignment_block <- c('BEGIN DATA;',
                       paste0('\tDIMENSIONS  NTAX=', n_sequences, ' NCHAR=',  n_positions, ';'),
                       '\tFORMAT DATATYPE=DNA GAP=- MISSING=?;',
                       '\tMATRIX',
                       labelled_alignment,
                       'END;',
                       '')
  
  
  # Format Assumptions box
  # This block inserts tip date precision 
  nexus_dates <- metadata_in_aln %>%
    select(tipnames, starts_with('date')) %>%
    
    # Set variable precision
    mutate(date_precision = case_when(
      date_format == 'yyyy' ~ 1,
      date_format == 'yyyy-mm' ~ 1/12,
      .default = 0)) %>%
    
    # for year-month and year only, date_dec = start
    mutate(date_lower = case_when(date_format == 'yyyy' ~ decimal_date(as_date(paste0(as.character(date_y), '-01-01'))),
                                  date_format == 'yyyy-mm' ~ decimal_date(as_date(paste0(date_ym, '-01'))),
                                  .default = date_dec)) %>%
    
    # sampling ranges
    mutate(date_upper = date_lower + date_precision) %>%
    
    # Write statements to be appended to nexus 
    mutate(nexus_statement = case_when(date_lower == date_upper ~ paste0('CALIBRATE ', tipnames, ' = fixed(', date_lower, '),'),
                                       .default = paste0('CALIBRATE ', tipnames, ' = uniform(', date_lower, ',', date_upper, '),'))) %>%
    pull(nexus_statement)
  
  # Write block
  assumptions_block <- list('BEGIN ASSUMPTIONS;',
                            '\tOPTIONS SCALE = years;',
                            paste0('\t' , dates[-length(dates)], collapse = '\n'),
                            paste0('\t' , gsub(',' , ';', dates[length(dates)]), collapse = '\n'),
                            'END;',
                            '')
  

  # If required, append the partitions for SRD06 model
  if(SRD06 == TRUE){
    # This block partitions the alignment by codon positions 1,2 and 3 as required for the 
    # SRD06 substitution model
    sets_block <- lc('BEGIN SETS;',
                     paste0('\tCHARSET 1,2 = 1-', ncol(aln)-2, '\\3, 2-', ncol(aln)-1, '\\3', ';'),
                     paste0('\tCHARSET 3 = 3-', ncol(aln), '\\3', ';'),
                     'END;',
                     '')
  }
  
  
  # Write file
  write_lines(filename, header_block)
  write_lines(filename, alignment_block, append = T)
  write_lines(filename, assumptions_block, append = T)
  
  if(SRD06 == TRUE){
    write_lines(filename, sets_block, append = T)
  }
  
  cat('file ', filename, ' written.')
}



WriteBEAST2Nexus(aln,
                 'test_2.nexus',
                 metadata = metadata)







lines_to_write <- c(header_block, 
                    data_block,
                    assumptions_block,
                    sets_block
                       )
write_lines(lines_to_write, 'test.txt')
names <- rownames(aln)

metadata_in_aln <- read_csv('./data/USUV_metadata_all_2025Feb10.csv') %>%
  filter(tipnames %in% names)
  



write_lines(dates, './2025Feb10/global_analysis/beast_mascot_constant/date_precision.txt')