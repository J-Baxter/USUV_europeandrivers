############# Dependencies ############# 
library(tidyverse)
library(ape)
library(ggmap)
library(sf)
library(giscoR)


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
data <-  read_csv('./data/ncbi_sequencemetadata_2024Aug13.csv', locale = readr::locale(encoding = "UTF-8")) %>%
  rows_update(updated_old_data %>%
                dplyr::select(Accession, Collection_Date, Host), 
              by = 'Accession', 
              unmatched = 'ignore') %>%
  left_join(updated_old_data %>% dplyr::select(-c(Collection_Date, Geo_Location, Host)), join_by(Accession)) %>% 
  dplyr::select(-c(ends_with('.y'), ends_with('.x')))  %>%
  # Allocate date format
  mutate(Collection_Date = gsub('[[:punct:]]', '-', Collection_Date)) %>%
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy",
    .default = 'missing')) 



# ANSES FRANCE Data
anses <- read_csv('./data/anses_data.csv', locale = readr::locale(encoding = "UTF-8")) %>%
  
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


nuts0 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "0")

nuts1 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "1")

nuts2 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "2")

nuts3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "3")

all_countries <- ne_countries(scale = 10, returnclass = "sf") %>%
  dplyr::select(name_en,
                iso_a2_eh,
                geometry) %>%
  rename(nuts0_id = iso_a2_eh,
         nuts0_name = name_en) %>%
  st_make_valid()

#############  Pipeline start ############# 
data_formatted_date <- data %>%
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
    .default = as.character(date_ymd)))

  
  # Format geolocation
data_formatted_geodata <- data_formatted_date %>%
  # Format available data
  mutate(across(starts_with('collection'), .fns = ~ tolower(.x))) %>%
  mutate(collection_location = str_trim(collection_location),
         collection_country = str_trim(collection_country)) %>%
  
  rowwise() %>%
  mutate(collection_tag = case_when(!grepl(collection_country, collection_location) ~ paste(collection_location, collection_country),
                                    .default = collection_location)) %>%
  #as_tibble() %>%
  
  # Using ggmaps:geocode, obtain precise lat-lon for each location. NB will give centroids
  # for large areas e.g cities or countries
  mutate(coords = geocode(collection_tag)) %>%
  as_tibble() %>%
  unnest(coords)
  
  # format as sf point
data_formatted_geodata_1 <- as_tibble(data_formatted_geodata) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  
  # Allocate NUTS codes for each sequence
  st_join(.,
          nuts0 %>%
            dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
          join = st_within, 
         # left = FALSE,
          largest = TRUE) %>%
  rename(geocode_coords = geometry,
         nuts0_id = NUTS_ID,
         nuts0_name = NUTS_NAME) %>%
  
  st_join(.,
          nuts1 %>%
            dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
          join = st_within, 
         # left = FALSE,
          largest = TRUE) %>%
  rename(nuts1_id = NUTS_ID,
         nuts1_name = NUTS_NAME) %>%
  
  st_join(.,
          nuts2 %>%
            dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
          join = st_within, 
          #left = FALSE,
          largest = TRUE) %>%
  rename(nuts2_id = NUTS_ID,
         nuts2_name = NUTS_NAME) %>%
  
  st_join(.,
          nuts3 %>%
            dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
          join = st_within, 
          #left = FALSE,
          largest = TRUE) %>%
  rename(nuts3_id = NUTS_ID,
         nuts3_name = NUTS_NAME) %>%
  
  
  # African and Asian Countries need some form of annnotation
  st_join(., 
          all_countries,
          join = st_within, 
          largest = TRUE) %>%
  mutate(nuts0_id = coalesce(nuts0_id.x, nuts0_id.y),
         nuts0_name = coalesce(nuts0_name.x, nuts0_name.y),
         .keep = 'unused')


  # Identify level of precision for each sequence: 
data_formatted_geodata_2 <- data_formatted_geodata_1 %>% 
  rowwise() %>%
  mutate(location_precision = case_when(
    #Country Only
    collection_country == collection_tag ~ 'nuts0',
    any(collection_location == tolower(nuts1$NUTS_NAME)) ~ 'nuts1',
    any(collection_location == tolower(nuts2$NUTS_NAME)) ~ 'nuts2',
    collection_country == 'greece' ~ 'nuts2',
    any(collection_location == tolower(nuts3$NUTS_NAME)) ~ 'nuts3',
    .default = 'exact')) %>%
  as_tibble() %>%
  base:: split(~location_precision) %>%
  map_at('nuts0', 
         ~ .x %>% mutate(across(starts_with("nuts1"), ~ NA),
                         across(starts_with("nuts2"), ~ NA),
                         across(starts_with("nuts3"), ~ NA))) %>%
  
  map_at('nuts1', 
         ~ .x %>% mutate(across(starts_with("nuts2"), ~ NA),
                         across(starts_with("nuts3"), ~ NA))) %>%
  
  map_at('nuts2', 
         ~ .x %>% mutate(across(starts_with("nuts3"), ~ NA))) %>%
  
  list_rbind()
  
  
data_formatted_host <- data_formatted_geodata_2 %>%
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
            by = join_by(primary_com_name)) 

# OU674388 - luxembourg date and host are in origial csv by now are NA?
  
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
  
data_formatted <- data_formatted_host %>%  
  # Arrange and rename columns
  dplyr::rename(
    #collection_regionname = region,
    #collection_countryname = country,
    #collection_countrycode = gid_0,
    #collection_countrylat = adm0_lat,
   #collection_countrylong = adm0_long,
    #collection_subdiv1name = name_1,
    #collection_subdiv1code = hasc_1,
    #collection_subdiv1lat = adm1_lat,
    #collection_subdiv1long = adm1_long,
    #collection_subdiv2name = name_2,
    #collection_subdiv2code = hasc_2,
    #collection_subdiv2lat = adm2_lat,
    #collection_subdiv2long = adm2_long,
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
    #collection_regionname,
    #collection_countryname,
    #collection_countrycode,
    #collection_countrylat,
    #collection_countrylong,
    #collection_subdiv1name,
    #collection_subdiv1code,
    #collection_subdiv1lat,
    #collection_subdiv1long,
    #collection_subdiv2name,
    #collection_subdiv2code,
    #collection_subdiv2lat,
    #collection_subdiv2long,
    nuts0_name,
    nuts0_id,
    nuts1_name,
    nuts1_id,
    nuts2_name,
    nuts2_id,
    nuts3_name,
    nuts3_id,
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
  dplyr::select(-c(organism_name,
            assembly,
            species,
            sequence_length,
            usa,
            collection_country,
            collection_date,
            notes,
            `migration pattern`,
            `wild/captive`,
            collection_tag,
            complete_location,
            complete_date,
            sequence_completeness,
            org_location,
            host,
            isolation_source)) %>%
  
  # location codes
  #mutate(collection_tipcode = coalesce(nuts3_id, nuts2_id, nuts1_id, nuts0_id)) %>%
  #mutate(unique_id = coalesce(sequence_accession, sequence_isolate)) %>%
  #mutate(collection_tipcode = gsub('\\.', '_', collection_tipcode)) %>%
  
  # Generate sequence names
  unite(.,
        tipnames, 
        sequence_accession, 
        sequence_isolate,
        host_sciname,
        nuts0_id,
        date_tipdate,
        sep = '|',                        
        remove = FALSE) %>%
  
  # remove non-permitted characters from tipnames
  mutate(tipnames = gsub(' ', '_', tipnames)) %>%
  mutate(tipnames = gsub('\\.', '', tipnames)) %>%
  mutate(tipnames = gsub('[^A-Za-z0-9_.\\|/\\-]', '_', tipnames)) %>%
  
  dplyr::select(-date_tipdate)



#############  Import alignment ############# 
alignment <- ape::read.dna('./2024Aug13/alignments/USUV_2024Aug13_alldata_aligned.fasta',
                           format = 'fasta',
                           as.matrix = T,
                           as.character = T) 


rownames(alignment) <- ReNameAlignment(alignment, data_formatted) 


start <- apply(alignment, 1, function(x) match(letters, x)[1])
end <- ncol(alignment) - apply(alignment[,ncol(alignment):1], 1, function(x) match(letters, x)[1])



sub_alignment_coords <- tibble(start = c(1003, 8968, 9100, 10042),
                               region = c('env', 'ns5', 'ns5', 'ns5'),
                               end = c(1491, 9264, 9600, 10314))


coords <- cbind('sequence_start' = start,
                           'sequence_end' = end) %>%
  as_tibble(rownames = 'tipnames') %>%
  mutate(sequence_length = sequence_end - sequence_start) %>%
  rowwise() %>%
  mutate(generegion_NS5_9000_9600 = ifelse(sequence_start <= 9100 && sequence_end >= 9150, 1, 0),
         #generegion_NS5_9100_9600 = ifelse(sequence_start <= 9200 && sequence_end >= 9500, 1, 0),
         generegion_NS5_10042_10312 = ifelse(sequence_start <= 10142 && sequence_end >= 10300, 1, 0),
         generegion_env_1003_1491 = ifelse(sequence_start <= 1100 && sequence_end >= 1200, 1, 0),
         generegion_nflg = ifelse(sequence_length >4200, 1, 0))

unnassigned <- coords %>%
  filter(if_all(starts_with('generegion'), ~ .x ==0))


#############  Previous Lineage assignments ############# 
old_lineages <- read_csv('./data/existing_lineages.csv')

#############  Join sequence data to metadata ############# 
metadata <- data_formatted %>%
  left_join(., coords) %>%
  left_join(old_lineages) %>%
  
  # Include only data in alignment
  #filter(tipnames %in% rownames(alignment)) %>%
  left_join(data %>% dplyr::select(Accession, Isolate, Release_Date), 
            by = join_by(sequence_accession == Accession,
                         sequence_isolate == Isolate)) %>%
  # Key to exclude FLI
  mutate(drop_fli = case_when(grepl('Friedrich-Loeffler-Institut', organization) & dmy(Release_Date) > as_date("2020-01-01") ~ TRUE,
                              .default = FALSE)) %>%
  dplyr::select(-Release_Date)  %>%
  
  #Europe
  mutate(is_europe = case_when(nuts0_id %in% nuts0$NUTS_ID ~ 1, .default = 0))


 
#############  Write metadata to file ############# #
write_csv(metadata, './data/USUV_metadata_all_2024Oct20.csv')


#############  Write alignment to file ############# 
alignment <- alignment[order(start),]
write.dna(alignment, './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted.fasta', format = 'fasta')


#############  Alignment without FLI ############# 
metadata_noFLI <- metadata %>%
  filter(!drop_fli)

write.dna(alignment[rownames(alignment) %in% metadata_noFLI$tipnames,], './2024Oct20/alignments/USUV_2024Oct20_alldata_aligned_formatted_noFLI.fasta', format = 'fasta')
write_csv(metadata_noFLI , './data/USUV_metadata_noFLI_2024Oct20.csv')


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
