# Initial Data Collection
library(rentrez)
library(tidyverse)
library(parsedate)
library(XML)
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
  #set_entrez_key(key)
  
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


MyFun3 <- function(country_code, region){
  
  if(!exists('iso2')){
    iso2 <- read_csv('./data/iso3166-2.csv') %>%
      select(-country_code) %>%
      separate_wider_delim(code, delim = '-', names = c('country_code', 'subdividsion_code'))
  }
  
  code <- iso2 %>%
    filter(country_code == country_code) %>%
    filter(subdivision_name == region) %>%
    pull(subdividsion_code) %>%
    unique()
  
  if(identical(code, character(0))){
    code <- NA
  }else
    
    return(list(code))
  
  
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
  mutate(country = countrycode::countrycode(country, 'country.name',
                                            destination = 'iso2c')) %>%
  mutate(region = case_when(
    
    #AT
    grepl('Burgenland|Stegersbach|Mattersburg|Bad Sauerbrunn', region) & country == 'AT' ~ 'Burgenland',
    
    grepl('Karnten', region) & country == 'AT' ~ 'Karnten',
    
    grepl('Niederosterreich|Voesendorf|Strasshof|Stockerau|Sollenau|Breitensee|Biberbach|Sankt Poelten|
          Pottschach|Neunkirchen|Markgrafneusiedl|Koenigstetten|Klosterneuburg|Hof am Leithaberge|
          Haringsee|Fishamend|Bruderndorf', region) & country == 'AT' ~ 'Niederosterreich',
    
    grepl('Oberosterreich|Sankt Georgen|Linz', region) & country == 'AT' ~ 'Oberosterreich',
    
    grepl('Salzburg', region) & country == 'AT' ~ 'Salzburg',
    
    grepl('Steiermark|Ragnitz|Graz|Eggersdorf', region) & country == 'AT' ~ 'Steiermark',
    
    grepl('Tirol', region) & country == 'AT' ~ 'Tirol',
    
    grepl('Vorarlberg', region) & country == 'AT' ~ 'Vorarlberg',
    
    grepl('Wien|Vienna', region) & country == 'AT' ~ 'Wien',
    
    
    #CZ
    grepl('Jihocesky kraj|Lomnice nad Luznici', region) & country == 'CZ' ~ 'Jihocesky kraj',
    
    grepl('Jihomoravsky kraj|Breclav|Brno', region) & country == 'CZ' ~ 'Jihomoravsky kraj',
    
    grepl('Karlovarsky kraj', region) & country == 'CZ' ~ 'Karlovarsky kraj',
    
    grepl('Kraj Vysocina', region) & country == 'CZ' ~ 'Kraj Vysocina',
    
    grepl('Kralovehradecky kraj', region) & country == 'CZ' ~ 'Kralovehradecky kraj',
    
    grepl('Liberecky kraj', region) & country == 'CZ' ~ 'Liberecky kraj',
    
    grepl('Moravskoslezsky kraj', region) & country == 'CZ' ~ 'Moravskoslezsky kraj',
    
    grepl('Olomoucky kraj', region) & country == 'CZ' ~ 'Olomoucky kraj',
    
    grepl('Pardubicky kraj', region) & country == 'CZ' ~ 'Pardubicky kraj',
    
    grepl('Plzensky kraj', region) & country == 'CZ' ~ 'Plzensky kraj',
    
    grepl('Praha, Hlavni mesto', region) & country == 'CZ' ~ 'Praha, Hlavni mesto',
    
    grepl('Stredocesky kraj|Prague', region) & country == 'CZ' ~ 'Stredocesky kraj',
    
    grepl('Ustecky kraj', region) & country == 'CZ' ~ 'Ustecky kraj',
    
    grepl('Zlinsky kraj', region) & country == 'CZ' ~ 'Zlinsky kraj',
    
    
    
    # IT
    grepl('Abruzzo', region) & country == 'IT' ~ 'Abruzzo',
    
    grepl('Basilicata', region) & country == 'IT'  ~ 'Basilicata',
    
    grepl('Basilicata', region) & country == 'IT'  ~ 'Calabria',
    
    grepl('Basilicata', region) & country == 'IT'  ~ 'Campania',
    
    grepl('Emilia|Bologna|Ferrara|Modena', region) & country == 'IT'  ~ 'Emilia-Romagna',
    
    grepl('Friuli Venezia Giulia|GO|PN|UD', region) ~ 'Friuli-Venezia Giulia',
    
    grepl('Lazio|Latium', region) & country == 'IT'  ~ 'Lazio',
    
    grepl('Liguria', region) & country == 'IT'  ~ 'Liguria',
    
    grepl('Lombard', region) & country == 'IT'  ~ 'Lombardia',
    
    grepl('Marche|Ancona|Macerata|Monteprandone-Ascoli Piceno|Osimo|Pesaro|Senigallia', region) & 
      country == 'IT'  ~ 'Marche',
    
    grepl('Molise', region) & country == 'IT'  ~ 'Molise',
    
    grepl('Piemonte|Casale Monferrato|Piedmont|Giarole|Mombello Monferrato|Verbania', region) & 
      country == 'IT'  ~ 'Piemonte',
    
    grepl('Puglia', region) & country == 'IT'  ~ 'Puglia',
    
    grepl('Sardegna', region) & country == 'IT'  ~ 'Sardegna',
    
    grepl('Sicilia', region) & country == 'IT'  ~ 'Sicilia',
    
    grepl('Toscana', region) & country == 'IT'  ~ 'Toscana',
    
    grepl('Trentino-Alto Adige|TN', region) & country == 'IT'  ~ 'Trentino-Alto Adige',
    
    grepl('Umbria', region) & country == 'IT'  ~ 'Umbria',
    
    grepl("Valle d'Aosta", region) & country == 'IT'  ~ "Valle d'Aosta",
    
    grepl('Veneto|VI|PD|\\bRO\\b|TV|VE|VR', region) & country == 'IT' ~ 'Veneto',
    
    
    # SN
    grepl('Barkedji', region) & country == 'SN' ~ 'Louga',
    
    
    #NL
    grepl('Drenthe|Westerveld|Borger-Odoorn|Tynaarlo|Eext|Koekange|MiddenDrexhe|Noordenveld', region) & 
      country == 'NL' ~ 'Drenthe',
    
    grepl('Flevoland|Almere|Lelystad|Noordoostpolder|Noordosterpolder', region) & 
      country == 'NL'~ 'Flevoland',
    
    grepl('Fryslan|Tytsjerksteradiel|Eastermar|Heerenveen', region) & country == 'NL'~ 'Fryslan',
    
    grepl('Gelderland|Zutphen|Zevenaar|Westvoort|Westervoort|Andelst|Arnhem|Bennekom|DeLiemers|
          Spijk|Doetichem|Doetinchem|Epe|Ermelo|Gelmonde|Ingen|Klarenbeek|Lochem|Oldenbroek', region) &
      country == 'NL' ~ 'Gelderland',
    
    grepl('Groningen|Zuidhorn|Oldambt|Spijk', region) & country == 'NL' ~ 'Groningen',
    
    grepl('Limburg|Venlo|Gennep|Kerkrade|Landgraaf|Ottersum', region) & country == 'NL' ~ 'Limburg',
    
    grepl('Noord-Brabant|Wernhout|Vlijmen|BeekEnDonk|Best|Boekel|Uden|Agatha|Schaik|Rosmalen|
          Reek|Lierop|Oss', region) ~ 'Noord-Brabant',
    
    grepl('Noord-Holland|Bloemendaal|Heerhugowaard|Hilversum|Huizen|Naarden|Reek', region) &
      country == 'NL' ~ 'Noord-Holland',
    
    grepl('Overijssel|Wierden|Enschede|Hardenberg|Heino|Losser|Overdinkel|Raalte', region) & 
      country == 'NL' ~ 'Overijssel',
    
    grepl('Utrecht|Bilthoven|Bosch en Duin|Bunnik|DeBilt|DeRondeVenen|Soest|Haarzuilens|Houten|
          Ijsselstein|Langbroek|Nieuwegein|Ijsselstein', region) & country == 'NL' ~ 'Utrecht',
    
    grepl('Zeeland|Terneuzen|Reimerswaal|Grenspad|Middelburg', region) & country == 'NL' ~ 'Zeeland',
    
    grepl('Zuid-Holland|Zoeterwoude|Westland|Den Haag|Rotterdam|Rijswijk|Leiden|Nederlek|HardinxveldGiessendam', region) & 
      country == 'NL' ~ 'Zuid-Holland',
    
    
    #RS
    grepl('Sabac', region) & country == 'RS' ~ 'Macvanski okrug',
    
    
    #DE
    grepl('Baden-Wuerttemberg|Freiburg im Breisgau|Rosengarten|Stutensee', region) & 
      country == 'DE' ~ 'Baden-Wuerttemberg',
    
    grepl('Bavaria|Munich|Nuremberg|Wuerzburg', region) & country == 'DE' ~ 'Bayern',
    
    grepl('Berlin', region) & country == 'DE' ~ 'Berlin',
    
    grepl('Brandenburg', region) & country == 'DE' ~ 'Brandenburg',
    
    grepl('Bremen', region)  & country == 'DE'~ 'Bremen',
    
    grepl('Hamburg', region) & country == 'DE' ~ 'Hamburg',
    
    grepl('Hessen|Freigericht|Giessen|Karben|Ranstadt', region) & country == 'DE' ~ 'Hessen',
    
    grepl('Mecklenburg-Vorpommern|Grevesmuehlen|Parchim|Ruegen|Warin', region) &
      country == 'DE' ~ 'Mecklenburg-Vorpommern',
    
    grepl('Lower Saxony|Brietlingen|Dinklage|Hannover|Lueneburg|Nordhorn|Osnabruck|Osterholz-Scharnbeck
          |Theene|Wingst', region) & country == 'DE' ~ 'Niedersachsen',
    
    grepl('Aachen|Brueggen', region) & country == 'DE' ~ 'Nordrhein-Westfalen',
    
    grepl('Rheinland-Pfalz', region) & country == 'DE' ~ 'Rheinland-Pfalz',
    
    grepl('Saarland|Puettlingen', region) & country == 'DE' ~ 'Saarland',
    
    grepl('Sachsen|Doberschuetz|Dresden|Gross Dueben|Leipzig', region) & country == 'DE' ~ 'Sachsen',
    
    grepl('Sachsen-Anhalt|Halle|Zeitz', region) & country == 'DE' ~ 'Sachsen-Anhalt',
    
    grepl('Aumuehle|Brickeln|Luebeck|Pinneberg|Prohn', region) & country == 'DE' ~ 'Schleswig-Holstein',
    
    grepl('Thuringen|Gera', region) & country == 'DE' ~ 'Thuringen',
    
    
    #HU
    grepl('Bacs-Kiskun|Kecskemet|Kiskunfelegyhaza|Kunadacs', region) & country == 'HU' ~ 'Bacs-Kiskun',
    grepl('Baranya', region) & country == 'HU' ~ 'Baranya',
    grepl('Bekes', region) & country == 'HU' ~ 'Bekes',
    grepl('Borsod-Abauj-Zemplen', region) & country == 'HU' ~ 'Borsod-Abauj-Zemplen',
    grepl('Budapest', region) & country == 'HU' ~ 'Budapest',
    grepl('Csongrad-Csanad', region) & country == 'HU' ~ 'Csongrad-Csanad',
    grepl('Fejer', region) & country == 'HU' ~ 'Fejer',
    grepl('Gyor-Moson-Sopron', region) & country == 'HU' ~ 'Gyor-Moson-Sopron',
    grepl('Hajdu-Bihar', region) & country == 'HU' ~ 'Hajdu-Bihar',
    grepl('Heves', region) & country == 'HU' ~ 'Heves',
    grepl('Jasz-Nagykun-Szolnok', region) & country == 'HU' ~ 'Jasz-Nagykun-Szolnok',
    grepl('Komarom-Esztergom', region) & country == 'HU' ~ 'Komarom-Esztergom',
    grepl('Nograd', region) & country == 'HU' ~ 'Nograd',
    grepl('\\bPest\\b|\\bVac\\b', region) & country == 'HU' ~ 'Pest',
    grepl('Somogy|Kaposvar', region) & country == 'HU' ~ 'Somogy',
    grepl('Szabolcs-Szatmar-Bereg|Fabianhaza', region) & country == 'HU' ~ 'Szabolcs-Szatmar-Bereg',
    grepl('\\bTolna\\b', region) & country == 'HU' ~ 'Tolna',
    grepl('\\bVas\\b', region) & country == 'HU' ~ 'Vas',
    grepl('Veszprem', region) & country == 'HU' ~ 'Veszprem',
    grepl('Zala|Lovaszi|Nagylengyel|Zalaegerszeg', region) & country == 'HU' ~ 'Zala',
    
    
    # Serbia
    grepl('Beograd', region) & country == 'SR' ~ 'Beograd',
    
    grepl('Borski okrug', region) & country == 'SR' ~ 'Borski okrug',
    
    grepl('Branicevski okrug', region) & country == 'SR' ~ 'Branicevski okrug',
    
    grepl('Jablanicki okrug', region) & country == 'SR' ~ 'Jablanicki okrug',
    
    grepl('Juznobacki okrug', region) & country == 'SR' ~ 'Juznobacki okrug',
    
    grepl('Juznobanatski okrug', region) & country == 'SR' ~ 'Juznobanatski okrug',
    
    grepl('Kolubarski okrug', region) & country == 'SR' ~ 'Kolubarski okrug',
    
    grepl('Kosovsko-Mitrovacki okrug', region) & country == 'SR' ~ 'Kosovsko-Mitrovacki okrug',
    
    grepl('Macvanski okrug', region) & country == 'SR' ~ 'Macvanski okrug',
    
    grepl('Moravicki okrug', region) & country == 'SR' ~ 'Moravicki okrug',
    
    grepl('Nisavski okrug', region) & country == 'SR' ~ 'Nisavski okrug',
    
    grepl('Pcinjski okrug', region) & country == 'SR' ~ 'Pcinjski okrug',
    
    grepl('Pecki okrug', region) & country == 'SR' ~ 'Pecki okrug',
    
    grepl('Pirotski okrug', region) & country == 'SR' ~ 'Pirotski okrug',
    
    grepl('Podunavski okrug', region) & country == 'SR' ~ 'Podunavski okrug',
    
    grepl('Pomoravski okrug', region) & country == 'SR' ~ 'Pomoravski okrug',
    
    grepl('Prizrenski okrug', region) & country == 'SR' ~ 'Prizrenski okrug',
    
    grepl('Rasinski okrug', region) & country == 'SR' ~ 'Rasinski okrug',
    
    grepl('Raski okrug', region) & country == 'SR' ~ 'Raski okrug',
    
    grepl('Severnobacki okrug', region) & country == 'SR' ~ 'Severnobacki okrug',
    
    grepl('Severnobanatski okrug', region) & country == 'SR' ~ 'Severnobanatski okrug',
    
    grepl('Srednjebanatski okrug', region) & country == 'SR' ~ 'Srednjebanatski okrug',
    
    grepl('Sremski okrug', region) & country == 'SR' ~ 'Sremski okrug',
    
    grepl('Sumadijski okrug', region) & country == 'SR' ~ 'Sumadijski okrug',
    
    grepl('Toplicki okrug', region) & country == 'SR' ~ 'Toplicki okrug',
    
    grepl('Zajecarski okrug', region) & country == 'SR' ~ 'Zajecarski okrug',
    
    grepl('Zapadnobacki okrug', region) & country == 'SR' ~ 'Zapadnobacki okrug',
    
    grepl('Zlatiborski okrug', region) & country == 'SR' ~ 'Zlatiborski okrug',
    
    
    #Croatia - HR
    [1] "Bjelovarsko-bilogorska zupanija" "Brodsko-posavska zupanija"       "Dubrovacko-neretvanska zupanija"
    [4] "Grad Zagreb"                     "Istarska zupanija"               "Karlovacka zupanija"            
    [7] "Koprivnicko-krizevacka zupanija" "Krapinsko-zagorska zupanija"     "Licko-senjska zupanija"         
    [10] "Medimurska zupanija"             "Osjecko-baranjska zupanija"      "Pozesko-slavonska zupanija"     
    [13] "Primorsko-goranska zupanija"     "Sibensko-kninska zupanija"       "Sisacko-moslavacka zupanija"    
    [16] "Splitsko-dalmatinska zupanija"   "Varazdinska zupanija"            "Viroviticko-podravska zupanija" 
    [19] "Vukovarsko-srijemska zupanija"   "Zadarska zupanija"               "Zagrebacka zupanija" 
    
    "Zagreb"       "Jastrebarsko" "Nasic" 
    "Catalonia, Viladecans" #Spain
    "Ticin" #CH
    
    # Spain
    [1] "Andalucia"                   "Aragon"                      "Asturias, Principado de"    
    [4] "Canarias"                    "Cantabria"                   "Castilla y Leon"            
    [7] "Castilla-La Mancha"          "Catalunya"                   "Ceuta"                      
    [10] "Extremadura"                 "Galicia"                     "Illes Balears"              
    [13] "La Rioja"                    "Madrid, Comunidad de"        "Melilla"                    
    [16] "Murcia, Region de"           "Navarra, Comunidad Foral de" "Pais Vasco"                 
    [19] "Valenciana, Comunidad"      
    
    #CH
    [1] "Aargau"                 "Appenzell Ausserrhoden" "Appenzell Innerrhoden"  "Basel-Landschaft"      
    [5] "Basel-Stadt"            "Bern"                   "Fribourg"               "Geneve"                
    [9] "Glarus"                 "Graubunden"             "Jura"                   "Luzern"                
    [13] "Neuchatel"              "Nidwalden"              "Obwalden"               "Sankt Gallen"          
    [17] "Schaffhausen"           "Schwyz"                 "Solothurn"              "Thurgau"               
    [21] "Ticino"                 "Uri"                    "Valais"                 "Vaud"                  
    [25] "Zug"                    "Zurich"
    
    # FR
    [1] "Auvergne-Rhone-Alpes"       "Bourgogne-Franche-Comte"    "Bretagne"                  
    [4] "Centre-Val de Loire"        "Corse"                      "Grand-Est"                 
    [7] "Hauts-de-France"            "Ile-de-France"              "Normandie"                 
    [10] "Nouvelle-Aquitaine"         "Occitanie"                  "Pays-de-la-Loire"          
    [13] "Provence-Alpes-Cote-d'Azur"
    
    # Mopping up
    is.na(region) ~ NA,
    grepl('Devexer', region) ~ NA, 
    .default = region
    
  )) %>%
  rowwise() %>%
  mutate(subdivision = MyFun3(country, region)) %>%
  
  
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
    host.class,
    host.order,
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