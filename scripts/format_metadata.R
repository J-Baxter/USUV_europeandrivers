############# Dependencies ############# 
library(tidyverse)
library(ape)


#############  Required user functions ############# 
ReNameAlignment <- function(alignment, data){
  names <- rownames(alignment)
  isolates <- rep(NA, length(names))
  
  published_isolates <- which(str_count(names, '\\|') >=2)
  izsve_isolates <- which(grepl('\\d{2}(VIR|pool)\\d{3,4}', names))
  anses_isolates <- which(grepl('\\/France\\/', names))
  apha_isolates <- which(grepl('^Blackbird', names))
  
  isolates[published_isolates] <- str_extract(names[published_isolates], "([^|]*)\\||([^,]*),")%>%
    str_replace_all("\\||,", "") %>%
    str_trim() %>%
    gsub('\\..*', '', .) 
  
  
  isolates[izsve_isolates] <- gsub("_[A-Za-z].+", "", names[izsve_isolates])%>%
    str_trim()
  
  isolates[anses_isolates] <- str_extract(names[anses_isolates], "USUV-[A-Za-z]{0,8}\\d+\\/France\\/\\d{4}") %>%
    str_trim()
  
  isolates[apha_isolates] <- names[apha_isolates] %>%
    str_trim()
  
  
  new_seqnames <- c()
  for (i in 1:length(isolates)){
    new_seqnames[i] <- data$tipnames[which(data$sequence_accession == isolates[i] | data$sequence_isolate == isolates[i])]
    
  }
  
  return(new_seqnames)
}

#############  Required source files ############# 
source('./scripts/FormatBirds.R')
source('./scripts/FormatMammals.R')
source('./scripts/FormatMosquitos.R')
source('./scripts/FormatTicks.R')
source('./scripts/FormatGeo.R')


#############  Required data ############# 

# Import previous search (that contains manual edits) and update with data from FLI
old_data <- read_csv('./data/ncbi_sequencemetadata_Apr12.csv') %>%
  dplyr::select(c(Accession, Collection_Date, Geo_Location, Host)) 

fli_updated <- read_csv('./data/fli_supplementary/fli_dates.csv') %>%
  mutate(across(ends_with('Date'), .fns = ~dmy(.x))) %>%
  mutate(Admission_Date = format(Admission_Date, '%Y-%m')) %>%
  mutate(Collection_Date = coalesce(as.character(Collection_Date), Admission_Date)) %>%
  dplyr::select(-Admission_Date) %>%
  filter(Accession %in% old_data$Accession)

updated_old_data <- old_data %>%
  rows_update(fli_updated, by = 'Accession')


# NCBI Genbank Data
data <-  read_csv('./data/ncbi_sequencemetadata_2024Aug13.csv') %>%
  left_join(updated_old_data, join_by(Accession)) %>%
  mutate(Collection_Date = coalesce(Collection_Date.y, Collection_Date.x),
         Geo_Location = coalesce(Geo_Location.y, Geo_Location.x),
         Host = coalesce(Host.x, Host.y)) %>%
  dplyr::select(-c(ends_with('.y'), ends_with('.x')))  %>%
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd", 
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date) & as.numeric(str_split_i('date', '-', 2)) <= 12 ~ "dd-mm-yyyy", 
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm", 
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy" )) 


# ANSES FRANCE Data
anses <- read_csv('./data/anses_data.csv') %>%
  
  # Clean date col
  mutate(Collection_Date = case_when(
    grepl('appro{0,1}x.', Collection_Date) ~ gsub(' appro{0,1}x.', '', Collection_Date) %>% gsub('^[^/]*/', '', .),
    .default = Collection_Date) %>%
      gsub('/', '-', .)) %>%
  
  # Allocate date format
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy",
    .default = 'missing'))%>%
  
  # Unpublished data only
  filter(is.na(Accession))


# APHA Data
apha <- read_csv('./data/apha_data.csv') %>%
  
  # Clean date col
  mutate(Collection_Date = case_when(
    grepl('appro{0,1}x.', Collection_Date) ~ gsub(' appro{0,1}x.', '', Collection_Date) %>% gsub('^[^/]*/', '', .),
    .default = Collection_Date) %>%
      gsub('/', '-', .)) %>%
  
  # Allocate date format
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy",
    .default = 'missing')) %>%
  
  # Unpublished data only
  filter(is.na(Accession))

# Izsve Data
izsve <- read_csv('./data/izsve_data.csv') %>%
  
  # Clean date col
  mutate(Collection_Date = case_when(
    grepl('appro{0,1}x.', Collection_Date) ~ gsub(' appro{0,1}x.', '', Collection_Date) %>% gsub('^[^/]*/', '', .),
    .default = Collection_Date) %>%
      gsub('/', '-', .)) %>%
  
  # Allocate date format
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy",
    .default = 'missing')) %>%
  
  # Unpublished data only
  filter(is.na(Accession))


geodata <- read_csv('data/updated_geodata.csv')
birds <- read_csv('bird_taxonomy.csv')
mammals <- read_csv('mammal_taxonomy.csv')
ticks <- read_csv('tick_taxonomy.csv') %>%
  mutate(sci_name = scientificName) %>%
  rename(primary_com_name = scientificName) %>%
  dplyr::select(-genus)
mosquitos <- read_csv('mosquito_taxonomy.csv') %>%
  mutate(sci_name = scientificName) %>%
  rename(primary_com_name = scientificName) %>%
  dplyr::select(-genus)

all_taxa <- bind_rows(birds, mammals, mosquitos, ticks)


#############  Pipeline start ############# 
data_formatted <- data %>%
  bind_rows(.,anses) %>%
  bind_rows(.,apha) %>%
  bind_rows(.,izsve) %>%

  # format colnames
  rename_with(., ~ tolower(.x)) %>%
  dplyr::select(-c(submitters, release_date)) %>%
  rename(sequence_length = length,
         sequence_accession = accession,
         sequence_isolate = isolate,
         sequence_completeness = nuc_completeness,
         collection_country = country,
         collection_location = geo_location) %>%
  
  mutate(collection_location = coalesce(collection_location, collection_country)) %>%

  # format date (create date_dec, date_ym date_ymd) %>%
  mutate(collection_date = str_trim(collection_date))%>%
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', collection_date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', collection_date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', collection_date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', collection_date) ~ "mm-yyyy",
    grepl('\\d{4}', collection_date) ~ "yyyy",
    .default = 'missing')) %>%
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
  dplyr::select(-date_parsed) %>%

  mutate(date_tipdate = case_when(
    is.na(date_ymd) & is.na(date_ym) ~ date_y,
    is.na(date_ymd) & !is.na(date_ym) ~ date_ym,
    .default = as.character(date_ymd))) %>%
  #select(-c(collection_date, date_format, complete_date)) %>%
 
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
    collection_country == 'poland' && !is.na(collection_location) ~ FormatPoland(collection_location), 
    collection_country == 'greece' && !is.na(collection_location) ~ FormatGreece(collection_location), 
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
  mutate(primary_com_name = FormatTicks(primary_com_name)) %>%
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
 # select(-c(collection_location, 
           # collection_country,
           # host, 
           # `migration pattern`,
            ##notes, 
           # match,
            #iso_1, 
           # `wild/captive`, 
           # sequence_length,
           # complete_location)) %>%
  
  # location codes
  mutate(collection_tipcode = coalesce(collection_subdiv2code, collection_subdiv1code, collection_countrycode)) %>%
  mutate(unique_id = coalesce(sequence_accession, sequence_isolate)) %>%
  mutate(collection_tipcode = gsub('\\.', '_', collection_tipcode)) %>%
  
  # Generate sequence names
  unite(.,
        tipnames, 
        sequence_accession, 
        sequence_isolate,
        host_sciname,
        collection_tipcode,
        date_tipdate,
        sep = '|',                        
        remove = FALSE) %>%
  
  # remove non-permitted characters from tipnames
  mutate(tipnames = gsub(' ', '_', tipnames)) %>%
  mutate(tipnames = gsub('\\.', '', tipnames)) %>%
  mutate(tipnames = gsub('[^A-Za-z0-9_.\\|/\\-]', '_', tipnames)) %>%
  
  dplyr::select(-unique_id)



#############  Import alignment ############# 
alignment <- ape::read.dna('./2024Aug13/alignments/USUV_2024Aug13_alldata_aligned.fasta',
                           format = 'fasta',
                           as.matrix = T,
                           as.character = T) 


rownames(alignment) <- ReNameAlignment(alignment , data_formatted) 


start <- apply(alignment, 1, function(x) match(letters, x)[1])
end <- ncol(alignment) - apply(alignment[,ncol(alignment):1], 1, function(x) match(letters, x)[1])



sub_alignment_coords <- tibble(start = c(1003, 8968, 9100, 10042),
                               region = c('env', 'ns5', 'ns5', 'ns5'),
                               end = c(1491, 9264, 9600, 10314))


coords <- cbind('sequence_start' = start,
                           'sequence_end' = end) %>%
  as_tibble(rownames = 'tipnames') %>%
  mutate(sequence_length = sequence_end - sequence_start) %>%
  mutate(sequence_generegion = case_when(
    sequence_start < 400 & sequence_end >9500 ~ 'nflg',
    sequence_start < 1100 & sequence_end <3000 ~ 'env_1003_1491',
    sequence_start > 7000 & sequence_end < 9300 ~ 'NS5_8968_9264',
    sequence_start > 7000 & sequence_end > 9400 ~ 'NS5_9100_9600',
    sequence_start > 9000 & sequence_end > 10000 ~ 'NS5_10042_10314'
  )) 


#############  Join sequence data to metadata ############# 
metadata <- data_formatted %>%
  select(-c(
    organism_name,
    assembly,
    species,
    collection_country,
    collection_location,
    host,
    collection_date,
    collection_location,
    date_format,
    isolation_source,
    complete_date,
    complete_location,
    `wild/captive`,
    `migration pattern`,
    notes,
    match,
    iso_1,
    collection_tipcode,
    sequence_length,
    sequence_completeness,
    usa,
    tissue_specimen_source
  )) %>%
  left_join(., coords) %>%
  
  # Include only data in alignment
  filter(tipnames %in% rownames(alignment)) %>%
  left_join(data %>% select(Accession, Isolate, Organization, Release_Date), 
            by = join_by(sequence_accession == Accession,
                         sequence_isolate == Isolate)) %>%
  # Key to exclude FLI
  mutate(drop_fli = case_when(grepl('Friedrich-Loeffler-Institut', Organization) & Release_Date > as_date("2020-01-01") ~ TRUE,
                              .default = FALSE)) %>%
  select(-c(Release_Date)) %>%
  rename(collection_organsiation = Organization)


 
#############  Write metadata to file ############# #
write_csv(metadata, './data/USUV_metadata_all_2024Aug28.csv')


#############  Write alignment to file ############# 
alignment <- alignment[order(start),]
write.dna(alignment, './2024Aug13/alignments/USUV_2024Aug28_alldata_aligned_formatted.fasta', format = 'fasta')


#############  Alignment without FLI ############# 
metadata_noFLI <- metadata %>%
  filter(!drop_fli)

write.dna(alignment[rownames(alignment) %in% metadata_noFLI$tipnames,], './2024Aug13/alignments/USUV_2024Aug28_alldata_aligned_formatted_noFLI.fasta', format = 'fasta')
write_csv(metadata_noFLI , './data/USUV_metadata_noFLI_2024Aug28.csv')
