#register_google('AIzaSyBkz-A4kEs2LKFrjNtbLyaQdLtH_gedMF8')

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
  
  
  rowwise() %>%
  mutate(collection_tag = case_when(!grepl(collection_country, collection_location) ~ paste(collection_location, collection_country),
                                    .default = collection_location)) %>%
  as_tibble() %>%
  
  mutate(coords = geocode(collection_tag)) %>%
  
  unnest(coords) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  
  #join country
  st_join(nuts0 %>%
            dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
          join = st_within, 
          largest = TRUE, 
          left = FALSE) %>%
  rename(geocode_coords = geometry,
         nuts0_id = NUTS_ID,
         nuts0_name = NUTS_NAME) %>%
  
  # Identify level of precision for each sequence: 
  mutate(location_precision = case_when(
    collection_country == collection_tag ~ 'nuts0',
    collection_country == collection_tag ~ 'nuts1',
    collection_country == collection_tag ~ 'nuts2',
    collection_country == collection_tag ~ 'nuts3',
    collection_country == collection_tag ~ 'exact',
    .default = NA_character_
    
  ))
  mutate(country_only = case_when(collection_country == collection_tag ~ '1',
                                  .default = '0')) %>%
  
  # Identify sequences where only NUTS2 is known
  mutate(best_resolution = case_when(grepl() ~ '1',
                                  .default = '0')) %>%
  
  # Allocate NUTS2 labels
  split(~country_only) %>% 
  map_at('0',  
         ~ st_join(.x, 
                   nuts1 %>% 
                     dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
                   join = st_within,
                   largest = TRUE, 
                   left = FALSE) %>%
           rename(nuts1_code= NUTS_ID,
                  nuts1_name= NUTS_NAME) %>%
           
           st_join(nuts2 %>% 
                     dplyr::select(geometry, NUTS_ID, NUTS_NAME), 
                   join = st_within,
                   largest = TRUE, 
                   left = FALSE) %>%
           rename(nuts2_code= NUTS_ID,
                  nuts2_name= NUTS_NAME)) %>%
  map_at('1',  
         ~ mutate(.x, 
                  nuts1_code = NA_character_,
                  nuts1_name = NA_character_,
                  nuts2_code = NA_character_,
                  nuts2_name = NA_character_)) %>%
  list_rbind() 
           
