############# Dependencies ############# 
library(tidyverse)
library(ape)


#############  Required user functions ############# 
ReNameAlignment <- function(alignment, data){
  isolates <- str_extract(rownames(alignment), "([^|]*)\\||([^,]*),")%>%
    str_replace_all("\\||,", "") %>%
    str_trim() %>%
    gsub('\\..*', '', .) %>% 
    gsub('Usutu virus isolate ','', .)
  
  new_seqnames <- c()
  
  
  for (i in 1:length(isolates)){
    newname <- vector()
    
    if(any(data$sequence_accession %in% isolates[[i]])){
      newname <- data$tipnames[data$sequence_accession %in% isolates[[i]]]
      
    }else if(any(data$sequence_isolate %in% isolates[[i]])){
      newname <- data$tipnames[data$sequence_isolate %in% isolates[[i]]]
    }
    
    if (length(newname)>0){
      new_seqnames[i] <-  newname
    }else{
      print(i)
      new_seqnames[i] <- isolates[[i]]
    }
    
  }
  
  rownames(alignment) <-  new_seqnames
  
  return(alignment)
}

#############  Required source files ############# 
source('./scripts/FormatBirds.R')
source('./scripts/FormatMammals.R')
source('./scripts/FormatMosquitos.R')
source('./scripts/FormatGeo.R')


#############  Required data ############# 
data <-  read_csv('./data/ncbi_sequencemetadata_Apr12.csv')
anses <- read_csv('./data/anses_data.csv') %>%
  mutate(Collection_Date = case_when(
    grepl('appro{0,1}x.', Collection_Date) ~ gsub(' appro{0,1}x.', '', Collection_Date) %>% gsub('^[^/]*/', '', .),
    .default = Collection_Date) %>%
      gsub('/', '-', .)) %>%
  mutate(date_format = case_when(
    str_count(Collection_Date, "-") == 2 ~ "dd-mm-yyyy",
    str_count(Collection_Date, "-") == 1 ~ "mm-yyyy"))

geodata <- read_csv('data/updated_geodata.csv')
birds <- read_csv('bird_taxonomy.csv')
mammals <- read_csv('mammal_taxonomy.csv')
mosquitos <- read_csv('mosquito_taxonomy.csv') %>%
  mutate(sci_name = scientificName) %>%
  rename(primary_com_name = scientificName) %>%
  select(-genus)

all_taxa <- bind_rows(birds, mammals, mosquitos)


#############  Pipeline start ############# 
data_formatted <- data%>%
  bind_rows(.,anses) %>%

  # format colnames
  rename_with(., ~ tolower(.x)) %>%
  select(-c(submitters, isolation_source, organization, release_date)) %>%
  rename(sequence_length = length,
         sequence_accession = accession,
         sequence_isolate = isolate,
         sequence_completeness = nuc_completeness,
         collection_country = country,
         collection_location = geo_location) %>%
  
  mutate(collection_location = coalesce(collection_location, collection_country)) %>%

  # format date (create date_dec, date_ym date_ymd) %>%
  mutate(collection_date = str_trim(collection_date))%>%
  split(~date_format) %>% 
  map_at("yyyy-mm-dd",  
         ~ mutate(.x,
                  date_parsed = ymd(collection_date),
                  date_ymd = format(date_parsed, '%Y-%m-%d'), 
                  date_dec = decimal_date(date_parsed),
                  date_ym = format(date_parsed, '%Y-%m'),
                  date_y = format(date_parsed, '%Y') )) %>% 
  map_at("dd-mm-yyyy",  
         ~ mutate(.x,
                  date_parsed = dmy(collection_date),
                  date_ymd = format(date_parsed, '%Y-%m-%d'), 
                  date_dec = decimal_date(date_parsed),
                  date_ym = format(date_parsed, '%Y-%m'),
                  date_y = format(date_parsed, '%Y') )) %>% 
  map_at("yyyy-mm",
         ~ mutate(.x,
                  date_parsed = ym(collection_date),
                  date_ymd = NA_character_, 
                  date_dec = NA,
                  date_ym = format(date_parsed,'%Y-%m') ,
                  date_y = format(date_parsed, '%Y') )) %>%
  map_at("mm-yyyy",
         ~ mutate(.x,
                  date_parsed = my(collection_date),
                  date_ymd = NA_character_, 
                  date_dec = NA,
                  date_ym = format(date_parsed,'%Y-%m') ,
                  date_y = format(date_parsed, '%Y') )) %>%
  map_at("yyyy",
         ~ mutate(.x,
                  date_parsed = as.POSIXlt.character(collection_date, format = '%Y'),
                  date_ymd = NA_character_, 
                  date_ym = NA_character_,
                  date_dec = NA,
                  date_y = format(date_parsed, '%Y'))) %>%

  list_rbind() %>%
  #select(-date_parsed) %>%

  mutate(date_tipdate = case_when(
    is.na(date_ymd) & is.na(date_ym) ~ date_y,
    is.na(date_ymd) & !is.na(date_ym) ~ date_ym,
    .default = as.character(date_ymd))) %>%
  select(-c(collection_date, date_format, complete_date)) %>%
 
   # Format geolocation
  mutate(across(starts_with('collection'), .fns = ~ tolower(.x))) %>%
  mutate(collection_location = str_trim(collection_location)) %>%
  rowwise() %>%
  mutate(match = case_when(
    collection_country == 'austria' && !is.na(collection_location) ~ FormatAustria(collection_location),
    collection_country == 'belgium' && !is.na(collection_location) ~ FormatBelgium(collection_location), 
    collection_country == 'croatia' && !is.na(collection_location) ~ FormatCroatia(collection_location),
    collection_country == 'czech republic' && !is.na(collection_location) ~ FormatCzechia(collection_location), 
    collection_country == 'france' && !is.na(collection_location) ~ FormatFrance(collection_location),
    collection_country == 'germany' && !is.na(collection_location) ~ FormatGermany(collection_location),
    collection_country == 'hungary' && !is.na(collection_location) ~ FormatHungary(collection_location), 
    collection_country == 'italy' && !is.na(collection_location) ~ FormatItaly(collection_location),
    collection_country == 'netherlands' && !is.na(collection_location) ~ FormatNetherlands(collection_location),
    collection_country == 'romania' && !is.na(collection_location) ~ FormatRomania(collection_location),
    collection_country == 'senegal' && !is.na(collection_location) ~ FormatSenegal(collection_location), 
    collection_country == 'serbia' && !is.na(collection_location) ~ FormatSerbia(collection_location),
    collection_country == 'slovakia' && !is.na(collection_location) ~ FormatSlovakia(collection_location),
    collection_country == 'spain' && !is.na(collection_location) ~ FormatSpain(collection_location),
    collection_country == 'switzerland' && !is.na(collection_location) ~ FormatSwitzerland(collection_location),
    collection_country == 'uganda' && !is.na(collection_location) ~ FormatUganda(collection_location), 
    collection_country == 'uk' && !is.na(collection_location) ~ FormatUK(collection_location), 
    .default = collection_country)) %>%
  as_tibble() %>%
  left_join(geodata, 
            by = join_by(match == match),
            na_matches =  "never") %>%

  
  # Format host 
  mutate(host = gsub("_", " ", tolower(host))) %>%
  mutate(host = str_trim(host)) %>%
  rowwise() %>%
  mutate(primary_com_name = FormatBird(host)) %>%
  mutate(primary_com_name = FormatMammal(primary_com_name)) %>%
  mutate(primary_com_name = FormatMosquito(primary_com_name)) %>%
  as_tibble() %>%
  left_join(all_taxa, 
            by = join_by(primary_com_name)) %>%
  
  # Additional Host info: local/migratory/captive
  # Bergmann et al R = resident species, P = partial migrants, S = short-distance migrants, L = long-distance migrants, Captive, NA
 # mutate(host_migration = case_when(
   # primary_com_name %in% c('european blackbird', 'blue tit') ~ 'resident',
  #  primary_com_name %in% c() ~ 'partial',
   # primary_com_name %in% c() ~ 'short-distance',
   # primary_com_name %in% c() ~ 'long-distance',
   # primary_com_name %in% c() ~ 'captive',
   # .default = NA_character_
  #)) %>%
  
  
  # Arrange and rename columns
  dplyr::rename(
    collection_regionname = region,
    collection_countryname = country,
    collection_countrycode = gid_0,
    collection_countrylat = adm0_lat,
    collection_countrylong = adm0_long,
    collection_subdiv1name = name_1,
    collection_subdiv1code = hasc_1,
    collection_subdiv1lat = adm1_lat,
    collection_subdiv1long = adm1_long,
    collection_subdiv2name = name_2,
    collection_subdiv2code = hasc_2,
    collection_subdiv2lat = adm2_lat,
    collection_subdiv2long = adm2_long,
    host_class = class,
    host_order = order,
    host_family = family,
    host_sciname = sci_name,
    host_commonname = primary_com_name
    )  %>%
  
  dplyr::relocate(
    sequence_accession,
    sequence_isolate,
    sequence_length,
    sequence_completeness,
    collection_regionname,
    collection_countryname,
    collection_countrycode,
    collection_countrylat,
    collection_countrylong,
    collection_subdiv1name,
    collection_subdiv1code,
    collection_subdiv1lat,
    collection_subdiv1long,
    collection_subdiv2name,
    collection_subdiv2code,
    collection_subdiv2lat,
    collection_subdiv2long,
    date_ymd,
    date_dec,
    date_ym,
    date_y,
    date_tipdate,
    host_class,
    host_order,
    host_family,
    host_sciname,
    host_commonname) %>%
  
  # Format NAs
  mutate(across(everything(), .fns = ~ gsub('^NA$', NA, .x))) %>%
  
  # Select columns
  select(-c(collection_location, 
            collection_country,
            host, 
            `migration pattern`,
            notes, 
            match,
            iso_1, 
            `wild/captive`, 
            sequence_length,
            complete_location)) %>%
  
  # location codes
  mutate(collection_tipcode = coalesce(collection_subdiv2code, collection_subdiv1code, collection_countrycode)) %>%
  mutate(unique_id = coalesce(sequence_accession, sequence_isolate)) %>%
  mutate(collection_tipcode = gsub('\\.', '_', collection_tipcode)) %>%
  
  # Generate sequence names
  unite(.,
        tipnames, 
        unique_id, 
        host_sciname,
        collection_tipcode,
        date_tipdate,
        sep = '|',                        
        remove = FALSE) %>%
  select(-unique_id) %>%
  mutate(tipnames = gsub(' ', '_', tipnames)) %>%
  mutate(tipnames = gsub('\\.', '', tipnames))
  



#############  Import alignment ############# 
alignment <- ape::read.dna('./2024Apr21/alignments/master_alignment.fasta',
                           format = 'fasta',
                           as.matrix = T,
                           as.character = T) %>%
  ReNameAlignment(. , data_formatted) 


start <- apply(alignment, 1, function(x) match(letters, x)[1])
end <- ncol(alignment) - apply(alignment[,ncol(alignment):1], 1, function(x) match(letters, x)[1])

coords <- cbind.data.frame('sequence_start' = start,
                           'sequence_end' = end) %>%
  as_tibble(rownames = 'tipnames') %>%
  mutate(sequence_length = sequence_end - sequence_start) %>%
  mutate(sequence_generegion = case_when(
    sequence_start < 400 & sequence_end >9500 ~ 'nflg',
    sequence_start < 1100 & sequence_end <3000 ~ 'E',
    sequence_start >7000 & sequence_start <8999~ 'NS5a',
    sequence_start >9000  ~ 'NS5b'
    
  )) 


#############  Join sequence data to metadata ############# 
metadata <- data_formatted %>%
  left_join(., coords) 


#############  Write metadata to file ############# 
write_csv(metadata, './data/metadata_2024Apr30.csv')


#############  Write alignment to file ############# 
write.dna(alignment, './2024Apr21/alignments/master_alignment_renamed_may29.fasta', format = 'fasta')


