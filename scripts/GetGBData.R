# Initial Data Collection
library(rentrez)
library(genbankr)
library(tidyverse)
source('./scripts/apikey.R')

# Function to extract citation, pubmed id, sequencing tech and accession from ID
ProcessEntrezXML <- function(xml){
  
  # Takes parsed XML
  if(any(class(xml) %in% c( "XMLInternalDocument" ,"XMLAbstractDocument"))){
    xml <- xmlToList(xml)
  }
  stopifnot(class(xml) == 'list')

  # Citation
  num_references <- length(xml[["GBSeq_references"]])
  journal <- c()
  pubmed <- c()
  
  for (i in 1:num_references){
    tmp <- NULL
    tmp <- try(xml[["GBSeq_references"]][[i]][["GBReference_journal"]])
    
    # Journal References
    if(!is.null(tmp)){
      journal[i] <- xml[["GBSeq_references"]][[i]][["GBReference_journal"]]
      
    }else if(is.null(tmp)){
      journal[i] <- NA
    }
    
    #Pubmed ID
    tmp <- NULL
    tmp <- try(xml[["GBSeq_references"]][[i]][["GBReference_pubmed"]])
    
    if(!is.null(tmp)){
      pubmed[i] <- xml[["GBSeq_references"]][[i]][["GBReference_pubmed"]] %>%
        as.numeric()
      
    }else if(is.null(tmp)){
      pubmed[i] <- NA
    }
    
  }
  
  # Sequencing Technology
  seqtech<- NULL
  tmp <- NULL
  tmp <- try(xml[["GBSeq_comment"]])
  
  if(!is.null(tmp)){
    seqtech <- xml[["GBSeq_comment"]] %>% 
      str_split(.,';') %>%
      sapply(., function(x) keep(x, ~ grepl('Sequencing', .x))) %>% 
      unlist()
    
    if(length(seqtech > 0)){
      seqtech <- seqtech %>%
        str_split_1(., '::') %>% 
        keep(., ~ grepl('sanger|roche|454|illumina|oxford|pac|nanopore', .x, ignore.case = T)) %>%
        str_trim()
    } else{
      seqtech <- NULL
    }
    
  }else if(is.null(tmp)){
    seqtech <- NA
  }
  
  # Seq length
  seqlength <- NULL
  tmp <- NULL
  tmp <- try(xml[["GBSeq_length"]])
  
  if(!is.null(tmp)){
    seqlength <- xml[["GBSeq_length"]] %>% 
      as.numeric()
    
  }else if(is.null(tmp)){
    seqlength <- NA
  }
  
  
  # Accession (To link back with other data)
  tmp <- NULL
  tmp <- try(xml[["GBSeq_locus"]])
  
  if(!is.null(tmp)){
    accn <- xml[["GBSeq_locus"]] 
    
  }else if(is.null(tmp)){
    accn <- NA
  }
  
  # Ensure empty elements are represented as NA
  if(identical(accn, character(0))){
    accn <- NA
  }
  
  if(identical(seqtech, character(0))){
    seqtech <- NA
  }
  
  if(identical(pubmed, character(0))){
    pubmed <- NA
  }
  
  if(identical(journal, character(0))){
    journal <- NA
  }
  
  if(identical(seqlength, character(0))){
    seqlength <- NA
  }
  
  out <- list('accession' = accn,
              'sequencing.technology' = seqtech,
              'pubmed.id' = pubmed, 
              'citation' = journal,
              'seqlength' = seqlength)
  return(out)
}


MyFunc2 <- function(dataframe){
  labels <- pull(dataframe, 16) %>% str_split(. ,'\\|')
  values <- pull(dataframe, 17) %>% str_split(. ,'\\|')
  accession <- select(dataframe, accession) 
  
  out <- mapply(function(x,y) rbind(x) %>% 
                  as_tibble() %>%
                  setNames(y), values, labels) %>% 
    bind_rows() %>%
    select(-strain) %>%
    cbind.data.frame(accession) 
  
  return(out)
}


# Wrapper to various Rentrez functions to generate dataframe containing:
# accession, seqname, genbank upload dates, molecule type, source, accession, sequencing tech,
# pubmed id and citation
GetEntrez <- function(searchterm, key = API_KEY, max = 10){
  set_entrez_key(key)
  
  # RENTREZ search
  test_search <- entrez_search(db = 'nuccore',
                               term =  searchterm,
                               retmax = max)
  test_id <- test_search$ids
  
  if(length(test_id) < 251){
    
    test_summ <- entrez_summary(db="nuccore", 
                                id=test_id) %>%
      lapply(., function(x) keep(x, is.character)) %>%
      lapply(., as_tibble) %>%
      bind_rows()%>%
      rename(accession = caption) %>%
      mutate(uid = as.numeric(uid)) 
    
    test_xml <- entrez_fetch(db="nuccore",
                             id=test_id, 
                             rettype="xml",
                             parsed=TRUE) %>% 
      xmlToList()
    
    
  }else{
    test_summ <- split(test_id,
                       ceiling(seq_along(test_id) / 200)) %>%
      lapply(.,  entrez_summary, db="nuccore") %>%
      list_flatten() %>%
      lapply(., function(x) keep(x, is.character)) %>%
      lapply(., as_tibble) %>%
      bind_rows()%>%
      rename(accession = caption) 
    
    
    test_xml <- split(test_id,
                      ceiling(seq_along(test_id) / 200)) %>%
      lapply(., entrez_fetch, db="nuccore",  rettype="xml",
             parsed=TRUE) %>%
      lapply(., xmlToList) %>%
      list_flatten()
    
  }
  
  xml_data <- lapply(test_xml, ProcessEntrezXML)  %>% 
    bind_rows() %>%
    group_by(accession, sequencing.technology) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = id, 
                values_from = c(pubmed.id, citation),
                names_vary = 'fastest')
  
  tmp <- left_join(test_summ, xml_data, by = join_by(accession)) 
  out <- MyFunc2(tmp) %>%
    right_join(tmp, #preserve entrez_data
               by = join_by(accession))
  return(out)
}


# Get GB data using GenbankR package
# contains: organism, mol_type, isolate, host, country, collection_date, accession
# Not always available due to differences in annotations

GetGenBankR <- function(accn, key = API_KEY){
  
  #set_entrez_key(key)
  
  require(genbankr)
  
  formatted_accn <- GBAccession(accn)
  #set_entrez_key(key)
  
  try_accn <-  try(readGenBank(formatted_accn, 
                               partial = TRUE),
                   silent = T)
  
  if(!class(try_accn) == "try-error"){
    entry <- readGenBank(formatted_accn, 
                         partial = TRUE)
    
    data <- entry@sources@elementMetadata@listData %>%  keep(is.character) %>% 
      as_tibble() %>%
      mutate(accession = accn)
    
  }else{
    data <- NULL
  }
  
  return(data)
  
}

###################################################################################################
# Function to extract Genbank collection dates from NCBI nucleotide database
# Requires function NCBI API key, otherwise downloads will exceed maximum number of calls

entrez_data <- GetEntrez("usutu[All Fields] AND viruses[filter]", 
                         max = 2000) 

formatted_data <- entrez_data %>%
  # Remove uninformative cols
  select(-c(term, 
            flags, 
            biomol,
            topology,
            subtype,
            subname,
            segsetsize, 
            sourcedb,
            tech, 
            completeness,
            projectid,
            genome, 
            contains('assembly'),
            gb_acronym,
            collected_by,
            clone,
            identified_by,
            extra,
            geneticcode,
            strand,
            biosample,
            accessionversion,
            lat_lon,
            uid)) %>%
  mutate(notes = coalesce(note, note...2, note...6)) %>%
  select(-c(note, note...2, note...6)) %>%
  mutate(lab.passaged = case_when(!is.na(lab_host) ~ 1,
                                  .default = 0)) %>%
  select(-lab_host) %>%
  
  # format year (NB requires manual editing due to additional dates present in isolate names)
  mutate(collection.date.formatted = case_when(
    nchar(collection_date) == 4 ~ ymd(collection_date, truncated = 2L) %>% parsedate::parse_date() %m+% months(6),
    .default = parsedate::parse_date(collection_date))) %>%
  mutate(collection.date.decimal = decimal_date(collection.date.formatted)) %>%
  select(-collection_date) %>%
  
  # host (species, genus)
  mutate(host.genusspecies = case_when(
    grepl('\\bAnopheles\\b', host) ~ 'Anopheles_NA',
    grepl('merula|blackbird', host) ~ 'Turdus_merula',
    grepl('carrion crow', host) ~  'Corvus_corone',
    grepl('\\bcrow\\b', host) ~  'Corvus_NA',
    grepl('\\bColumba sp.\\b', host) ~ 'Columba_NA',
    grepl('Crocidura|shrew', host) ~ 'Crocidura_NA' ,
    grepl('\\bCulex\\b|Culex sp.|\\bCulex cf.\\b|Culex pipiens/torrentium pool', host) ~ 'Culex_NA',
    grepl('modestus', host) ~ 'Culex_modestus' ,
    grepl('pipiens', host) ~  'Culex_pipiens',
    grepl('Blue tit', host) ~ 'Cyanistes_caeruleus',
    grepl('eurasian jay|glandarius', host) ~  'Garrulus_glandarius',
    grepl('goshawk', host) ~  'Accipiter_gentilis',
    grepl('great tit', host) ~  'Parus_major',
    grepl('hooded crow', host) ~  'Corvus_cornix',
    grepl('Mastomys natalensis', host) ~ 'Mastomys_natalensis',
    grepl('Melanitta nigra', host) ~ 'Melanitta_nigra',
    grepl('Turdus philomelos', host) ~ 'Turdus_philomelos',
    grepl('Turdus pilaris', host) ~ 'Turdus_pilaris',
    grepl('tawny owl', host) ~ 'Strix_aluco',
    grepl('Sitta europaea', host) ~ 'Sitta_europaea',
    grepl('Panurus biarmicus', host) ~ 'Panurus_biarmicus',
    grepl('Parus caeruleus', host) ~ 'Parus_caeruleus',
    grepl('Passer domesticus', host) ~ 'Passer_domesticus',
    grepl('Passer montanus', host) ~ 'Passer_montanus ',
    grepl('Rattus rattus', host) ~ 'Rattus_rattus',
    grepl('spectacled warbler', host) ~ 'Sylvia_conspicillata',
    grepl('sparrowhawk', host) ~ 'Accipiter_nisus',
    grepl('Strix nebulosa', host) ~ 'Strix_nebulosa',
    grepl('vulgaris', host) ~ 'Sturnus_vulgaris',
    grepl('Sylviayatricapilla', host) ~ 'Sylvia_atricapilla',
    grepl('Sylviayatricapilla', host) ~ 'Sylvia_atricapilla',
    grepl('Sylviayatricapilla', host) ~ 'Sylvia_atricapilla',
    grepl('Tachyeres', host) ~ 'Tachyeres_NA',
    grepl('Serinus canaria', host) ~ 'Serinus_canaria',
    grepl('free-ranging wild birds|goose|Avian|gull|Ixodida|magpie|mosquito|Laridae|Turdidae|\\bowl\\b|swift|swallow', host) ~  'NA_NA',
    .default = gsub(' ', '_', host))) %>%
  separate_wider_delim(host.genusspecies, delim = '_', names = c('host.genus', 'host.species')) %>%
  
  # host order
  mutate(host.order = case_when(
    grepl('Anopheles|Culex|Aedes', host.genus) ~ 'Diptera',
    grepl('Homo', host.genus) ~ 'Primates',
    grepl('Gallus|Alectoris|Phasianus', host.genus) ~ 'Galliformes',
    grepl('Corvus|Passer|Turdus|Sturnus|Sylvia|Alauda|Cyanistes|Erithacus|Erythrura|Fringilla|Garrulus|Gracula|
          Leucopsar|Muscicapa|Panurus|Parus|Pica|Pyrrhula|Serinus|Sitta', host.genus) ~ 'Passeriformes',
    grepl('Accipiter', host.genus) ~ 'Accipitriformes',
    grepl('Strix|Asio|Surnia', host.genus) ~ 'Strigiformes',
    grepl('Alopochen|Anas|Branta|Tachyeres', host.genus) ~ 'Anseriformes',
    grepl('Caprimulgus', host.genus) ~ 'Caprimulgiformes',
    grepl('Ciconia', host.genus) ~ 'Ciconiiformes',
    grepl('Columba|Streptopelia', host.genus) ~ 'Columbiformes',
    grepl('Larus', host.genus) ~ 'Charadriiformes',
    grepl('Mastomys|Rattus', host.genus) ~ 'Rodentia',
    grepl('Pipistrellus', host.genus) ~ 'Chiroptera',
    grepl('Ixodes', host) ~ 'Parasitiformes',
    grepl('swift|swallow', host) ~ 'Passeriformes',
    grepl('goose', host) ~ 'Anseriformes', 
    grepl('gull', host) ~ 'Charadriiformes' )) %>%

  # host class 
  mutate(host.class = case_when(
    grepl('Gallifomes|Passeriformes|Accipitriformes|Strigiformes|Anseriformes|Columbiformes|
          Charadriiformes', host.order) ~ 'Aves',
    grepl('Rodentia|Primates|Chiroptera', host.order) ~ 'Mammalia',
    grepl('Diptera', host.order) ~ 'Insecta',
    grepl('Parasitiformes', host.order) ~ 'Arachnida')) %>%
  
  mutate(across(contains('host'), .fns = ~ case_when(grepl('\\bNA\\b', .x) ~ NA,
                                                     .default = .x))) %>%
  select(-host) %>%
  
  # Location
  # Region detail (regions need grouping/harmonising)
  separate_wider_delim(country, delim = ':', names =  c('country', 'region'), too_few = 'align_start') %>%
  mutate(across(c(country, region), .fns = ~ str_trim(.x))) %>%
  
  # section of genome analysed
  mutate(full.length = case_when(seqlength > 1000 ~ 'nflg',
                                 .default = 'partial')) %>%
  rename(isolation.source = isolation_source) %>%
  rename(seq.length = seqlength) %>%
  relocate(
    accession,
    organism,
    isolate,
    collection.date.formatted,
    collection.date.decimal,
    lab.passaged,
    country,
    region,
    genotype,
    host.order,
    host.class,
    host.genus,
    host.species, 
    isolation.source, 
    seq.length,
    full.length,
    sequencing.technology,
    title,
    pubmed.id_1,
    citation_1) 

filtered_data <- formatted_data %>%
  # only USUV 
  filter(organism == 'Usutu virus') %>%
  
  # Non passaged
  filter(lab.passaged == 0) %>%
  select(-c(strain, organism, lab.passaged, createdate, updatedate))

write.csv(filtered_data, './data/genbank_usuv_20231206.csv')

 





###################################################################################################