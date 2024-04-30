data <-  read_csv('./data/ncbi_sequencemetadata_Apr12.csv')
geodata <- read_csv('data/updated_geodata.csv')
birds <- read_csv('bird_taxonomy.csv')
mammals <- read_csv('mammal_taxonomy.csv')
mosquitos <- read_csv('mosquito_taxonomy.csv') %>%
  mutate(sci_name = scientificName) %>%
  rename(primary_com_name = scientificName) %>%
  select(-genus)

all_taxa <- bind_rows(birds, mammals, mosquitos)

data_formatted <- data %>%
  
  # format colnames
  rename_with(., ~ tolower(.x)) %>%
  select(-c(submitters, isolation_source, organization, release_date)) %>%
  rename(sequence_length = length,
         sequence_accession = accession,
         sequence_isolate = isolate,
         sequence_completeness = nuc_completeness,
         collection_country = country,
         collection_location = geo_location) %>%
  
  # format date (create date_dec, date_ym date_ymd) %>%
  
  mutate(date_ymd = case_when(complete_date == TRUE ~ suppressWarnings(parsedate::parse_date(collection_date)) %>% 
                                as.Date() )) %>%
  mutate(date_dec = case_when(complete_date == TRUE ~ suppressWarnings(format(round(decimal_date(date_ymd), 2), 
                               nsmall = 2)))) %>%
  mutate(date_ym = case_when(complete_date == TRUE | date_format == 'yyyy-mm'~ 
                               suppressWarnings(parsedate::parse_date(collection_date)) %>% 
                               format(., "%Y-%m"))) %>%
  mutate(date_y = suppressWarnings(format(parsedate::parse_date(collection_date), "%Y"))) %>%
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
    collection_country == 'Czech Republic' && !is.na(collection_location) ~ FormatCzechia(collection_location), 
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
  select(-c(collection_location, collection_country, host, `migration pattern`, notes, match, iso_1, `wild/captive`)) %>%
  
  # location codes
  mutate(collection_tipcode = coalesce(collection_subdiv2code, collection_subdiv1code, collection_countrycode)) %>%
  mutate(collection_tipcode = gsub('\\.', '_', collection_tipcode)) %>%
  
  # Generate sequence names
  unite(.,
        tipnames, 
        sequence_accession, 
        host_sciname,
        collection_tipcode,
        date_tipdate,
        sep = '|',                        
        remove = FALSE) %>%
  mutate(tipnames = gsub(' ', '_', tipnames)) %>%
  mutate(tipnames = gsub('\\.', '', tipnames))
  

write_csv(data_formatted, './data/metadata_2024Apr30.csv')
