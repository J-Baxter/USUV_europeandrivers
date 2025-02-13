####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 202
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)
library(ggmap)
library(sf)
library(giscoR)
library(rnaturalearth)

# User functions
ReNameAlignment <- function(alignment, data){
  names <- rownames(alignment)
  isolates <- rep(NA, length(names))
  
  published_isolates <- which(str_count(names, '\\|') >=2)
  izsve_isolates <- which(grepl('^\\d{2}(VIR|pool)\\d{3,4}', names))
  anses_isolates <- which(grepl('^Usutu virus isolate', names))
  apha_isolates <- which(grepl('^Blackbird', names))
  erasmus_isolates <- which(grepl('^[a-zA-Z0-9_.-]{8,18}\\|\\d{4}(-\\d{2}-\\d{2}){0,1}', names))
  greece_isolates <- which(grepl('USUV-GR1439', names))
  
  isolates[published_isolates] <- str_extract(names[published_isolates], "([^|]*)\\||([^,]*),")%>%
    str_replace_all("\\||,", "") %>%
    str_trim() %>%
    gsub('\\..*', '', .) 
  
  
  isolates[izsve_isolates] <- gsub("_[A-Za-z].+", "", names[izsve_isolates])%>%
    str_trim()
  
  isolates[anses_isolates] <- str_extract(names[anses_isolates], "\\d+\\/France\\/\\d{2}\\/\\d{4}") %>%
    str_trim()
  
  isolates[apha_isolates] <- names[apha_isolates] %>%
    str_trim()
  
  isolates[erasmus_isolates] <- names[erasmus_isolates] %>%
    str_replace_all("\\|.*", "") %>%
    str_trim()
  
  isolates[greece_isolates] <- names[greece_isolates] %>%
    str_trim()
  
  
  new_seqnames <- c()
  for (i in 1:length(isolates)){

    new_name <- data$tipnames[which(data$sequence_accession == isolates[i] | data$sequence_isolate == isolates[i])]
    old_name <- names[i]
    
    if(length(new_name) > 0){
      new_seqnames[i] <- new_name
      
    }else{
      new_seqnames[i] <- old_name
    }
    
  }
  
  return(new_seqnames)
}


FormatNewData <- function(data){
  out <- data %>%
    
    # Clean date col
    mutate(Collection_Date = as.character(Collection_Date)) %>%
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
  
  return(out)
}


# Source Files
source('./scripts/FormatBirds.R')
source('./scripts/FormatMammals.R')
source('./scripts/FormatMosquitos.R')
source('./scripts/FormatTicks.R')


############################################## DATA ################################################
old_data <- read_csv('./data/ncbi_sequencemetadata_Apr12.csv') %>%
  dplyr::select(c(Accession, Collection_Date, Geo_Location, Host)) 

# FLI sequences with additional location/date information
fli_updates <- read_csv('./data/fli_supplementary/fli_dates.csv') %>%
  mutate(across(ends_with('Date'), .fns = ~dmy(.x))) %>%
  mutate(Admission_Date = format(Admission_Date, '%Y-%m')) %>%
  mutate(Collection_Date = coalesce(as.character(Collection_Date), Admission_Date)) %>%
  dplyr::select(-Admission_Date) %>%
  filter(Accession %in% old_data$Accession)

# Additional data from Daniel Cadar
cadar_metadata <- read_csv('./data/Dcadar_metadata.csv') %>%
  dplyr::select(accession, geo_location, collection_date, host) %>%
  rename_with(~gsub('_', ' ', .x) %>%
                str_to_title() %>%
                gsub(' ', '_', .)) %>%
  mutate(Collection_Date = dmy(Collection_Date) %>%
           format(., '%Y-%m-%d'))


# Import search from NCBI 2024-08-13
ncbi_aug13 <- read_csv('./data/ncbi_sequencemetadata_2024Aug13.csv', 
                       locale = readr::locale(encoding = "UTF-8")) 


# Import search from NCBI 2024-09 - 2024-12
ncbi_dec  <- read_csv('./ncbi_aug24todec24/ncbi_aug24-dec24_metadata.csv', 
                      locale = readr::locale(encoding = "UTF-8")) %>%
  mutate(Release_Date = as.character(Release_Date))

# Import search from NCBI 2024-09 - 2024-12
ncbi_feb  <- read_csv('./ncbi_dec24tofeb25/ncbi_dec24-feb25_metadata.csv', 
                      locale = readr::locale(encoding = "UTF-8")) %>%
  mutate(Release_Date = as.character(Release_Date))

# ANSES FRANCE Data
anses <- read_csv('./data/anses_data_updated.csv',
                  locale = readr::locale(encoding = "UTF-8")) %>%
  FormatNewData(.)


# APHA Data
apha <- read_csv('./data/apha_data.csv', 
                 locale = readr::locale(encoding = "UTF-8")) %>%
  FormatNewData(.)

# Izsve Data
izsve <- read_csv('./data/izsve_data.csv', 
                  locale = readr::locale(encoding = "UTF-8")) %>%
  FormatNewData(.)

# ERASMUS Data
erasmus <- read_csv('./data/erasmus_data.csv', 
                    locale = readr::locale(encoding = "UTF-8")) %>%
  FormatNewData(.)

# Greece Data
greece <- read_csv('./greece_data/greece_data.csv', 
                   locale = readr::locale(encoding = "UTF-8")) %>%
  FormatNewData(.)


# Portugal Data
#portugal <- read_csv('./portugal_data/portugal_data.csv', 
                  #   locale = readr::locale(encoding = "UTF-8")) %>%
  #FormatNewData(.)


# Import reference datasets for animal taxa
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


# Import shapefiles for geographic regions
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


############################################## MAIN ################################################

# Update previously maintained data with new FLI data
updated_old_data <- old_data %>%
  rows_update(fli_updates, by = 'Accession') 

# Merge 13th Aug data with old data
data <- ncbi_aug13 %>%
  rows_update(updated_old_data %>%
                dplyr::select(Accession, Collection_Date, Host), 
              by = 'Accession', 
              unmatched = 'ignore') %>%
  left_join(updated_old_data %>% dplyr::select(-c(Collection_Date, Geo_Location, Host)), join_by(Accession)) %>% 
  dplyr::select(-c(ends_with('.y'), ends_with('.x')))  %>%
  
  # Update with Cadar
  rows_update(cadar_metadata, by = 'Accession') %>%
  
  # Add New data post Dec 2024
  rows_insert(ncbi_dec) %>%
  rows_insert(ncbi_feb %>% mutate(Collection_Date = as.character(Collection_Date))) %>%
  
  # Allocate date format
  mutate(Collection_Date = gsub('[[:punct:]]', '-', Collection_Date)) %>%
  mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', Collection_Date)  ~ "yyyy-mm-dd",
    grepl('\\d{2}-\\d{2}-\\d{4}', Collection_Date)  ~ "dd-mm-yyyy",
    grepl('\\d{4}-\\d{2}', Collection_Date) ~ "yyyy-mm",
    grepl('\\d{2}-\\d{4}', Collection_Date) ~ "mm-yyyy",
    grepl('\\d{4}', Collection_Date) ~ "yyyy",
    .default = 'missing'))  


# <------- Main pipeline starts here
data_formatted_date <- data %>%
  bind_rows(.,anses) %>%
  bind_rows(.,apha) %>%
  bind_rows(.,izsve) %>%
  bind_rows(., erasmus) %>%
  bind_rows(., greece) %>%
  #bind_rows(., portugal) %>%
  
  
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
  
  # drop duplicate ZA sequences
  filter(!grepl('MF374485|EU074021|AY453412', sequence_accession)) %>%
  
  # drop reference sequence
  filter(!grepl('NC_006551', sequence_accession)) %>%
  
  # drop lab sequences (mostly senegal and france)
  filter(!grepl('PP482814|PP482817|LP944231|EU303236|LC227579|KU760915', sequence_accession)) %>%
  
  # drop KC754958 (highly divergent 1969 sequence)
  filter(!grepl('KC754958', sequence_accession)) %>%
  
  
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

data_formatted <- data_formatted_host %>%  
  # Arrange and rename columns
  dplyr::rename(
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




# Import alignment 
alignment <- ape::read.dna('./2025Feb10/alignments/USUV_2025Feb10_alldata_aligned.fasta',
                           format = 'fasta',
                           as.matrix = T,
                           as.character = T) 


rownames(alignment) <- ReNameAlignment(alignment, data_formatted) 

alignment <- alignment[rownames(alignment) %in% data_formatted$tipnames,]
agtc <- c('a', 't', 'g', 'c', 'x')
start <- apply(alignment, 1, function(x) match(agtc, x)[1])
end <- ncol(alignment) - apply(alignment[,ncol(alignment):1], 1, function(x) match(agtc, x)[1])

n_ambig <- apply(alignment, 1, function(x) str_c(x, collapse = '') %>%
                   sub(".*?[atcg]", "", .) %>%
                   str_count(., "n")) 


sub_alignment_coords <- tibble(start = c(1003, 8968, 9100, 10042),
                               region = c('env', 'ns5', 'ns5', 'ns5'),
                               end = c(1491, 9264, 9600, 10314))


coords <- cbind('sequence_start' = start,
                'sequence_end' = end) %>%
  as_tibble(rownames = 'tipnames') %>%
  mutate(sequence_length = sequence_end - sequence_start) %>%
  
  # ambiguities
  mutate(sequence_ambig = n_ambig) %>%
  mutate(sequence_ambig = n_ambig/sequence_length) %>%
  rowwise() %>%
  mutate(generegion_NS5 = ifelse(sequence_start <= 9100 && sequence_end >= 10300, 1, 0),
         generegion_env = ifelse(sequence_start <= 1100 && sequence_end >= 1200, 1, 0),
         generegion_nflg = ifelse(sequence_length > 9500, 1, 0))

unnassigned <- coords %>%
  filter(if_all(starts_with('generegion'), ~ .x ==0))


# Previous Lineage assignments #
old_lineages <- read_csv('./data/existing_lineages.csv')

##  Join sequence data to metadata # 
metadata <- data_formatted %>%
  left_join(., coords) %>%
  left_join(old_lineages) %>%
  
  # Include only data in alignment
  #filter(tipnames %in% rownames(alignment)) %>%
  left_join(data %>% dplyr::select(Accession, Isolate, Release_Date), 
            by = join_by(sequence_accession == Accession,
                         sequence_isolate == Isolate)) %>%
  # Key to exclude FLI
  mutate(drop_fli = case_when(grepl('Friedrich-Loeffler-Institut', organization) & dmy(Release_Date)  > as_date("2020-01-01") ~ TRUE,
                              .default = FALSE)) %>% #ignore warnings - these are because some Release_Date are NA or not in the correct format
  dplyr::select(-Release_Date)  %>%
  
  #Europe
  mutate(is_europe = case_when(nuts0_id %in% nuts0$NUTS_ID ~ 1, .default = 0))


############################################## WRITE ###############################################

# Write metadata to file #
write_csv(metadata, './data/USUV_metadata_all_2025Feb10.csv')


# Write alignment to file # 
alignment <-  alignment %>%
  as.DNAbin() %>%
  .[order(start),]

write.FASTA(alignment, './2025Feb10/alignments/USUV_2025Feb10_alldata_aligned_formatted.fasta')


# Alignment without FLI # 
metadata_noFLI <- metadata %>%
  filter(!drop_fli)

write.FASTA(alignment[rownames(alignment) %in% metadata_noFLI$tipnames,], './2025Feb10/alignments/USUV_2025Feb10_alldata_aligned_formatted_noFLI.fasta')


# NFLG alignments only # 
write.FASTA(alignment[rownames(alignment) %in% (metadata %>% 
                                                  filter(generegion_nflg == 1) %>% 
                                                  pull(tipnames)),],
            './2025Feb10/alignments/USUV_2025Feb10_alldata_aligned_formatted_NFLG.fasta')

write.FASTA(alignment[rownames(alignment) %in% (metadata_noFLI %>% 
                                                filter(generegion_nflg == 1) %>% 
                                                pull(tipnames)),],
          './2025Feb10/alignments/USUV_2025Feb10_alldata_aligned_formatted_noFLI_NFLG.fasta')


############################################## END #################################################
####################################################################################################
####################################################################################################