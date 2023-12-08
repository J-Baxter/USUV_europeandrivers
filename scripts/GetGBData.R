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



###################################################################################################
# Function to extract Genbank collection dates from NCBI nucleotide database
# Requires function NCBI API key, otherwise downloads will exceed maximum number of calls

entrez_data <- GetEntrez("usutu[All Fields] AND viruses[filter]", 
                         max = 2000) 

# import ISO data
iso2 <- read_csv('./data/iso3166-2.csv') %>%
  select(-country_code) %>%
  separate_wider_delim(code, delim = '-', names = c('iso.country.code', 'iso.subdivision.code')) %>%
  rename(iso.subdivision.name = subdivision_name)



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
    nchar(collection_date) == 4 ~ ymd(collection_date, truncated = 2L) %>%
      parsedate::parse_date() , #%m+% months(6)
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
    
    grepl('Culex sp|\\bCulex cf\\b|Culex pipiens/torrentium pool', host, perl = T) ~ 'Culex_NA',
    
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
    
    grepl('Tachyeres', host) ~ 'Tachyeres_NA',
    
    grepl('Serinus canaria', host) ~ 'Serinus_canaria',
    
    grepl('free-ranging wild birds|goose|Avian|gull|Ixodida|magpie|mosquito|Laridae|Turdidae|\\bowl\\b|swift|swallow', host) ~  'NA_NA',
    
    .default = gsub(' ', '_', host))) %>%
  separate_wider_delim(host.genusspecies, delim = '_', names = c('host.genus', 'host.species'), too_few = 'align_start') %>%
  
  # host order
  mutate(host.order = case_when(
    grepl('Anopheles|Culex|Aedes', host.genus) ~ 'Diptera',
    
    grepl('Homo', host.genus) ~ 'Primates',
    
    grepl('Gallus|Alectoris|Phasianus', host.genus) ~ 'Galliformes',
    
    grepl('Corvus|Passer|Turdus|Sturnus|Sylvia|Alauda|Cyanistes|Erithacus|Erythrura|Fringilla|Garrulus|Gracula|Leucopsar|Muscicapa|Panurus|Parus|Pica|Pyrrhula|Serinus|Sitta',
          host.genus) ~ 'Passeriformes',
    
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
    grepl('Gallifomes|Passeriformes|Accipitriformes|Strigiformes|Anseriformes|Columbiformes|Charadriiformes', host.order) ~ 'Aves',
    
    grepl('Rodentia|Primates|Chiroptera', host.order) ~ 'Mammalia',
    
    grepl('Diptera', host.order) ~ 'Insecta',
    
    grepl('Parasitiformes', host.order) ~ 'Arachnida')) %>%
  
  mutate(across(contains('host'), .fns = ~ case_when(grepl('\\bNA\\b', .x) ~ NA,
                                                     .default = .x))) %>%
  select(-host) %>%
  
  # Location
  # Region detail (regions need grouping/harmonising)
  separate_wider_delim(country, delim = ':', names =  c('iso.country', 'region'), too_few = 'align_start') %>%
  mutate(across(c(iso.country, region), .fns = ~ str_trim(.x))) %>%
  mutate(iso.country.code = countrycode::countrycode(iso.country, 'country.name',
                                            destination = 'iso2c')) %>%
  mutate(iso.subdivision.name = case_when(
    
    #AUSTRIA -A T
    grepl('Burgenland|Stegersbach|Mattersburg|Bad Sauerbrunn', region) & iso.country.code == 'AT' ~ 'Burgenland',
    
    grepl('Karnten', region) & iso.country.code == 'AT' ~ 'Karnten',
    
    grepl('Niederosterreich|Voesendorf|Strasshof|Stockerau|Sollenau|Breitensee|Biberbach|Sankt Poelten|Pottschach|Neunkirchen|Markgrafneusiedl|Koenigstetten|Klosterneuburg|Hof am Leithaberge|Haringsee|Fishamend|Bruderndorf',
          region) & iso.country.code == 'AT' ~ 'Niederosterreich',
    
    grepl('Oberosterreich|Sankt Georgen|Linz', region) & iso.country.code == 'AT' ~ 'Oberosterreich',
    
    grepl('Salzburg', region) & iso.country.code == 'AT' ~ 'Salzburg',
    
    grepl('Steiermark|Ragnitz|Graz|Eggersdorf', region) & iso.country.code == 'AT' ~ 'Steiermark',
    
    grepl('Tirol', region) & iso.country.code == 'AT' ~ 'Tirol',
    
    grepl('Vorarlberg', region) & iso.country.code == 'AT' ~ 'Vorarlberg',
    
    grepl('Wien|Vienna', region) & iso.country.code == 'AT' ~ 'Wien',
    
    
    #CZECHIA - CZ
    grepl('Jihocesky kraj|Lomnice nad Luznici', region) & iso.country.code == 'CZ' ~ 'Jihocesky kraj',
    
    grepl('Jihomoravsky kraj|Breclav|Brno|Hlohovec', region) & iso.country.code == 'CZ' ~ 'Jihomoravsky kraj',
    
    grepl('Karlovarsky kraj', region) & iso.country.code == 'CZ' ~ 'Karlovarsky kraj',
    
    grepl('Kraj Vysocina', region) & iso.country.code == 'CZ' ~ 'Kraj Vysocina',
    
    grepl('Kralovehradecky kraj', region) & iso.country.code == 'CZ' ~ 'Kralovehradecky kraj',
    
    grepl('Liberecky kraj', region) & iso.country.code == 'CZ' ~ 'Liberecky kraj',
    
    grepl('Moravskoslezsky kraj', region) & iso.country.code == 'CZ' ~ 'Moravskoslezsky kraj',
    
    grepl('Olomoucky kraj', region) & iso.country.code == 'CZ' ~ 'Olomoucky kraj',
    
    grepl('Pardubicky kraj', region) & iso.country.code == 'CZ' ~ 'Pardubicky kraj',
    
    grepl('Plzensky kraj', region) & iso.country.code == 'CZ' ~ 'Plzensky kraj',
    
    grepl('Praha, Hlavni mesto', region) & iso.country.code == 'CZ' ~ 'Praha, Hlavni mesto',
    
    grepl('Stredocesky kraj|Prague', region) & iso.country.code == 'CZ' ~ 'Stredocesky kraj',
    
    grepl('Ustecky kraj', region) & iso.country.code == 'CZ' ~ 'Ustecky kraj',
    
    grepl('Zlinsky kraj', region) & iso.country.code == 'CZ' ~ 'Zlinsky kraj',
    
    
    
    # ITALY - IT
    grepl('Abruzzo', region) & iso.country.code == 'IT' ~ 'Abruzzo',
    
    grepl('Basilicata', region) & iso.country.code == 'IT'  ~ 'Basilicata',
    
    grepl('Basilicata', region) & iso.country.code == 'IT'  ~ 'Calabria',
    
    grepl('Basilicata', region) & iso.country.code == 'IT'  ~ 'Campania',
    
    grepl('Emilia|Bologna|Ferrara|Modena|EmiliaRomagna', region) & iso.country.code == 'IT'  ~ 'Emilia-Romagna',
    
    grepl('Friuli Venezia Giulia|GO|PN|UD', region) ~ 'Friuli-Venezia Giulia',
    
    grepl('Lazio|Latium', region) & iso.country.code == 'IT'  ~ 'Lazio',
    
    grepl('Liguria', region) & iso.country.code == 'IT'  ~ 'Liguria',
    
    grepl('Lombard', region) & iso.country.code == 'IT'  ~ 'Lombardia',
    
    grepl('Marche|Ancona|Macerata|Monteprandone-Ascoli Piceno|Osimo|Pesaro|Senigallia', region) & 
      iso.country.code == 'IT'  ~ 'Marche',
    
    grepl('Molise', region) & iso.country.code == 'IT'  ~ 'Molise',
    
    grepl('Piemonte|Casale Monferrato|Piedmont|Giarole|Mombello Monferrato|Verbania', region) & 
      iso.country.code == 'IT'  ~ 'Piemonte',
    
    grepl('Puglia', region) & iso.country.code == 'IT'  ~ 'Puglia',
    
    grepl('Sardegna', region) & iso.country.code == 'IT'  ~ 'Sardegna',
    
    grepl('Sicilia', region) & iso.country.code == 'IT'  ~ 'Sicilia',
    
    grepl('Toscana', region) & iso.country.code == 'IT'  ~ 'Toscana',
    
    grepl('Trentino-Alto Adige|TN', region) & iso.country.code == 'IT'  ~ 'Trentino-Alto Adige',
    
    grepl('Umbria', region) & iso.country.code == 'IT'  ~ 'Umbria',
    
    grepl("Valle d'Aosta", region) & iso.country.code == 'IT'  ~ "Valle d'Aosta",
    
    grepl('Veneto|VI|PD|\\bRO\\b|TV|VE|VR', region) & iso.country.code == 'IT' ~ 'Veneto',
    
    
    # SENEGAL - SN
    grepl('Barkedji', region) & iso.country.code == 'SN' ~ 'Louga',
    
    
    # NETHERLANDS - NL
    grepl('Drenthe|Westerveld|Borger-Odoorn|Tynaarlo|Eext|Koekange|MiddenDrexhe|Noordenveld', region) & 
      iso.country.code == 'NL' ~ 'Drenthe',
    
    grepl('Flevoland|Almere|Lelystad|Noordoostpolder|Noordosterpolder', region) & 
      iso.country.code == 'NL'~ 'Flevoland',
    
    grepl('Fryslan|Tytsjerksteradiel|Eastermar|Heerenveen', region) & iso.country.code == 'NL'~ 'Fryslan',
    
    grepl('Gelderland|Zutphen|Zevenaar|Westvoort|Westervoort|Andelst|Arnhem|Bennekom|DeLiemers|Spijk|Doetichem|Doetinchem|Epe|Ermelo|Gelmonde|Ingen|Klarenbeek|Lochem|Oldenbroek',
          region) &
      iso.country.code == 'NL' ~ 'Gelderland',
    
    grepl('Groningen|Zuidhorn|Oldambt|Spijk', region) & iso.country.code == 'NL' ~ 'Groningen',
    
    grepl('Limburg|Venlo|Gennep|Kerkrade|Landgraaf|Ottersum', region) & iso.country.code == 'NL' ~ 'Limburg',
    
    grepl('Noord-Brabant|Wernhout|Vlijmen|BeekEnDonk|Best|Boekel|Uden|Agatha|Schaik|Rosmalen|Reek|Lierop|Oss',
          region) ~ 'Noord-Brabant',
    
    grepl('Noord-Holland|Bloemendaal|Heerhugowaard|Hilversum|Huizen|Naarden|Reek', region) &
      iso.country.code == 'NL' ~ 'Noord-Holland',
    
    grepl('Overijssel|Wierden|Enschede|Hardenberg|Heino|Losser|Overdinkel|Raalte', region) & 
      iso.country.code == 'NL' ~ 'Overijssel',
    
    grepl('Utrecht|Bilthoven|Bosch en Duin|Bunnik|DeBilt|DeRondeVenen|Soest|Haarzuilens|Houten|Ijsselstein|Langbroek|Nieuwegein|Ijsselstein', 
          region) & iso.country.code == 'NL' ~ 'Utrecht',
    
    grepl('Zeeland|Terneuzen|Reimerswaal|Grenspad|Middelburg', region) & iso.country.code == 'NL' ~ 'Zeeland',
    
    grepl('Zuid-Holland|Zoeterwoude|Westland|Den Haag|Rotterdam|Rijswijk|Leiden|Nederlek|HardinxveldGiessendam', region) & 
      iso.country.code == 'NL' ~ 'Zuid-Holland',
 
    
    # GERMANY -DE
    grepl('Baden-Wuerttemberg|Freiburg im Breisgau|Rosengarten|Stutensee', region) & 
      iso.country.code == 'DE' ~ 'Baden-Wuerttemberg',
    
    grepl('Bavaria|Munich|Nuremberg|Wuerzburg', region) & iso.country.code == 'DE' ~ 'Bayern',
    
    grepl('Berlin', region) & iso.country.code == 'DE' ~ 'Berlin',
    
    grepl('Brandenburg', region) & iso.country.code == 'DE' ~ 'Brandenburg',
    
    grepl('Bremen', region)  & iso.country.code == 'DE'~ 'Bremen',
    
    grepl('Hamburg', region) & iso.country.code == 'DE' ~ 'Hamburg',
    
    grepl('Hessen|Freigericht|Giessen|Karben|Ranstadt', region) & iso.country.code == 'DE' ~ 'Hessen',
    
    grepl('Mecklenburg-Vorpommern|Grevesmuehlen|Parchim|Ruegen|Warin', region) &
      iso.country.code == 'DE' ~ 'Mecklenburg-Vorpommern',
    
    grepl('Lower Saxony|Brietlingen|Dinklage|Hannover|Lueneburg|Nordhorn|Osnabruck|Osterholz-Scharnbeck|Theene|Wingst',
          region) & iso.country.code == 'DE' ~ 'Niedersachsen',
    
    grepl('Aachen|Brueggen', region) & iso.country.code == 'DE' ~ 'Nordrhein-Westfalen',
    
    grepl('Rheinland-Pfalz', region) & iso.country.code == 'DE' ~ 'Rheinland-Pfalz',
    
    grepl('Saarland|Puettlingen', region) & iso.country.code == 'DE' ~ 'Saarland',
    
    grepl('Sachsen|Doberschuetz|Dresden|Gross Dueben|Leipzig', region) & iso.country.code == 'DE' ~ 'Sachsen',
    
    grepl('Sachsen-Anhalt|Halle|Zeitz', region) & iso.country.code == 'DE' ~ 'Sachsen-Anhalt',
    
    grepl('Aumuehle|Brickeln|Luebeck|Pinneberg|Prohn', region) & iso.country.code == 'DE' ~ 'Schleswig-Holstein',
    
    grepl('Thuringen|Gera', region) & iso.country.code == 'DE' ~ 'Thuringen',
    
    
    # HUNGARY-HU
    grepl('Bacs-Kiskun|Kecskemet|Kiskunfelegyhaza|Kunadacs', region) & iso.country.code == 'HU' ~ 'Bacs-Kiskun',
    
    grepl('Baranya', region) & iso.country.code == 'HU' ~ 'Baranya',
    
    grepl('Bekes', region) & iso.country.code == 'HU' ~ 'Bekes',
    
    grepl('Borsod-Abauj-Zemplen', region) & iso.country.code == 'HU' ~ 'Borsod-Abauj-Zemplen',
    
    grepl('Budapest', region) & iso.country.code == 'HU' ~ 'Budapest',
    
    grepl('Csongrad-Csanad', region) & iso.country.code == 'HU' ~ 'Csongrad-Csanad',
    
    grepl('Fejer', region) & iso.country.code == 'HU' ~ 'Fejer',
    
    grepl('Gyor-Moson-Sopron|\\bGyor\\b', region) & iso.country.code == 'HU' ~ 'Gyor-Moson-Sopron',
    
    grepl('Hajdu-Bihar', region) & iso.country.code == 'HU' ~ 'Hajdu-Bihar',
    
    grepl('Heves', region) & iso.country.code == 'HU' ~ 'Heves',
    
    grepl('Jasz-Nagykun-Szolnok', region) & iso.country.code == 'HU' ~ 'Jasz-Nagykun-Szolnok',
    
    grepl('Komarom-Esztergom', region) & iso.country.code == 'HU' ~ 'Komarom-Esztergom',
    
    grepl('Nograd', region) & iso.country.code == 'HU' ~ 'Nograd',
    
    grepl('\\bPest\\b|\\bVac\\b', region) & iso.country.code == 'HU' ~ 'Pest',
    
    grepl('Somogy|Kaposvar', region) & iso.country.code == 'HU' ~ 'Somogy',
    
    grepl('Szabolcs-Szatmar-Bereg|Fabianhaza', region) & iso.country.code == 'HU' ~ 'Szabolcs-Szatmar-Bereg',
    
    grepl('\\bTolna\\b', region) & iso.country.code == 'HU' ~ 'Tolna',
    
    grepl('\\bVas\\b', region) & iso.country.code == 'HU' ~ 'Vas',
    
    grepl('Veszprem', region) & iso.country.code == 'HU' ~ 'Veszprem',
    
    grepl('Zala|Lovaszi|Nagylengyel|Zalaegerszeg', region) & iso.country.code == 'HU' ~ 'Zala',
    
    
    # SERBIA-RS
    grepl('Beograd', region) & iso.country.code == 'RS' ~ 'Beograd',
    
    grepl('Borski okrug', region) & iso.country.code == 'RS' ~ 'Borski okrug',
    
    grepl('Branicevski okrug', region) & iso.country.code == 'RS' ~ 'Branicevski okrug',
    
    grepl('Jablanicki okrug', region) & iso.country.code == 'RS' ~ 'Jablanicki okrug',
    
    grepl('Juznobacki okrug', region) & iso.country.code == 'RS' ~ 'Juznobacki okrug',
    
    grepl('Juznobanatski okrug', region) & iso.country.code == 'RS' ~ 'Juznobanatski okrug',
    
    grepl('Kolubarski okrug', region) & iso.country.code == 'RS' ~ 'Kolubarski okrug',
    
    grepl('Kosovsko-Mitrovacki okrug', region) & iso.country.code == 'RS' ~ 'Kosovsko-Mitrovacki okrug',
    
    grepl('Macvanski okrug|Sabac', region) & iso.country.code == 'RS' ~ 'Macvanski okrug',
    
    grepl('Moravicki okrug', region) & iso.country.code == 'RS' ~ 'Moravicki okrug',
    
    grepl('Nisavski okrug', region) & iso.country.code == 'RS' ~ 'Nisavski okrug',
    
    grepl('Pcinjski okrug', region) & iso.country.code == 'RS' ~ 'Pcinjski okrug',
    
    grepl('Pecki okrug', region) & iso.country.code == 'RS' ~ 'Pecki okrug',
    
    grepl('Pirotski okrug', region) & iso.country.code == 'RS' ~ 'Pirotski okrug',
    
    grepl('Podunavski okrug', region) & iso.country.code == 'RS' ~ 'Podunavski okrug',
    
    grepl('Pomoravski okrug', region) & iso.country.code == 'RS' ~ 'Pomoravski okrug',
    
    grepl('Prizrenski okrug', region) & iso.country.code == 'RS' ~ 'Prizrenski okrug',
    
    grepl('Rasinski okrug', region) & iso.country.code == 'RS' ~ 'Rasinski okrug',
    
    grepl('Raski okrug', region) & iso.country.code == 'RS' ~ 'Raski okrug',
    
    grepl('Severnobacki okrug', region) & iso.country.code == 'RS' ~ 'Severnobacki okrug',
    
    grepl('Severnobanatski okrug', region) & iso.country.code == 'RS' ~ 'Severnobanatski okrug',
    
    grepl('Srednjebanatski okrug', region) & iso.country.code == 'RS' ~ 'Srednjebanatski okrug',
    
    grepl('Sremski okrug', region) & iso.country.code == 'RS' ~ 'Sremski okrug',
    
    grepl('Sumadijski okrug', region) & iso.country.code == 'RS' ~ 'Sumadijski okrug',
    
    grepl('Toplicki okrug', region) & iso.country.code == 'RS' ~ 'Toplicki okrug',
    
    grepl('Zajecarski okrug', region) & iso.country.code == 'RS' ~ 'Zajecarski okrug',
    
    grepl('Zapadnobacki okrug', region) & iso.country.code == 'RS' ~ 'Zapadnobacki okrug',
    
    grepl('Zlatiborski okrug', region) & iso.country.code == 'RS' ~ 'Zlatiborski okrug',
    
    
    #Croatia - HR
    grepl('Bjelovarsko-bilogorska zupanija', region) & iso.country.code == 'HR' ~ 'Bjelovarsko-bilogorska zupanija',
    
    grepl('Brodsko-posavska zupanija', region) & iso.country.code == 'HR' ~ 'Brodsko-posavska zupanija',
    
    grepl('Grad Zagreb|\\bZagreb\\b', region) & iso.country.code == 'HR' ~ 'Grad Zagreb',
    
    grepl('Istarska zupanija', region) & iso.country.code == 'HR' ~ 'Istarska zupanija',
    
    grepl('Karlovacka zupanija', region) & iso.country.code == 'HR' ~ 'Karlovacka zupanija',
    
    grepl('Koprivnicko-krizevacka zupanija', region) & iso.country.code == 'HR' ~ 'Koprivnicko-krizevacka zupanija',
    
    grepl('Krapinsko-zagorska zupanija', region) & iso.country.code == 'HR' ~ 'Krapinsko-zagorska zupanija',
    
    grepl('Licko-senjska zupanija zupanija', region) & iso.country.code == 'HR' ~ 'Licko-senjska zupanijazupanija',
    
    grepl('Medimurska zupanija', region) & iso.country.code == 'HR' ~ 'Medimurska zupanija',
    
    grepl('Osjecko-baranjska zupanija', region) & iso.country.code == 'HR' ~ 'Osjecko-baranjska zupanija',
    
    grepl('Pozesko-slavonska zupanija', region) & iso.country.code == 'HR' ~ 'Pozesko-slavonska zupanija',
    
    grepl('Primorsko-goranska zupanija', region) & iso.country.code == 'HR' ~ 'Primorsko-goranska zupanija',
    
    grepl('Sibensko-kninska zupanija', region) & iso.country.code == 'HR' ~ 'Sibensko-kninska zupanija',
    
    grepl('Sisacko-moslavacka zupanija', region) & iso.country.code == 'HR' ~ 'Sisacko-moslavacka zupanija',
    
    grepl('Splitsko-dalmatinska zupanija', region) & iso.country.code == 'HR' ~ 'Splitsko-dalmatinska zupanija',
    
    grepl('Varazdinska zupanija', region) & iso.country.code == 'HR' ~ 'Varazdinska zupanija',
    
    grepl('Viroviticko-podravska zupanija|Nasic', region) & iso.country.code == 'HR' ~ 'Viroviticko-podravska zupanija',
    
    grepl('Vukovarsko-srijemska zupanija', region) & iso.country.code == 'HR' ~ 'Vukovarsko-srijemska zupanija',
    
    grepl('Zadarska zupanija', region) & iso.country.code == 'HR' ~ 'Zadarska zupanija',
    
    grepl('Zagrebacka zupanija|Jastrebarsko', region) & iso.country.code == 'HR' ~ 'Zagrebacka zupanija',
    
    
    # SPAIN - ES
    grepl('Andalucia', region) & iso.country.code == 'ES' ~ 'Andalucia',
    
    grepl('Aragon', region) & iso.country.code == 'ES' ~ 'Aragon',
    
    grepl('Asturias, Principado de', region) & iso.country.code == 'ES' ~ 'Asturias, Principado de',
    
    grepl('Canarias', region) & iso.country.code == 'ES' ~ 'Canarias',
    
    grepl('Cantabria', region) & iso.country.code == 'ES' ~ 'Cantabria',
    
    grepl('Castilla y Leon', region) & iso.country.code == 'ES' ~ 'Castilla y Leon',
    
    grepl('Castilla-La Mancha', region) & iso.country.code == 'ES' ~ 'Castilla-La Mancha',
    
    grepl('Catalunya|Catalonia', region) & iso.country.code == 'ES' ~ 'Catalunya',
    
    grepl('Ceuta', region) & iso.country.code == 'ES' ~ 'Ceuta',
    
    grepl('Extremadura', region) & iso.country.code == 'ES' ~ 'Extremadura',
    
    grepl('Galicia', region) & iso.country.code == 'ES' ~ 'Galicia',
    
    grepl('Illes Balears', region) & iso.country.code == 'ES' ~ 'Illes Balears',
    
    grepl('La Rioja', region) & iso.country.code == 'ES' ~ 'La Rioja',
    
    grepl('Madrid, Comunidad de', region) & iso.country.code == 'ES' ~ 'Madrid, Comunidad de',
    
    grepl('Melilla', region) & iso.country.code == 'ES' ~ 'Melilla',
    
    grepl('Murcia, Region de', region) & iso.country.code == 'ES' ~ 'Murcia, Region de',
    
    grepl('Navarra, Comunidad Foral de', region) & iso.country.code == 'ES' ~ 'Navarra, Comunidad Foral de',
    
    grepl('Pais Vasco', region) & iso.country.code == 'ES' ~ 'Pais Vasco',
    
    grepl('Valenciana, Comunidad', region) & iso.country.code == 'ES' ~ 'Valenciana, Comunidad',

    
    # Switzerland - CH
    grepl('Aargau', region) & iso.country.code == 'CH' ~ 'Aargau',
    
    grepl('Appenzell Ausserrhoden', region) & iso.country.code == 'CH' ~ 'Appenzell Ausserrhoden',
    
    grepl('Appenzell Innerrhoden', region) & iso.country.code == 'CH' ~ 'Appenzell Innerrhoden',
    
    grepl('Basel-Landschaft', region) & iso.country.code == 'CH' ~ 'Basel-Landschaft',
    
    grepl('Basel-Stadt', region) & iso.country.code == 'CH' ~ 'Basel-Stadt',
    
    grepl('Bern', region) & iso.country.code == 'CH' ~ 'Bern',
    
    grepl('Fribourg', region) & iso.country.code == 'CH' ~ 'Fribourg',
    
    grepl('Geneve', region) & iso.country.code == 'CH' ~ 'Geneve',
    
    grepl('Glarus', region) & iso.country.code == 'CH' ~ 'Glarus',
    
    grepl('Graubunden', region) & iso.country.code == 'CH' ~ 'Graubunden',
    
    grepl('Jura', region) & iso.country.code == 'CH' ~ 'Jura',
    
    grepl('Luzern', region) & iso.country.code == 'CH' ~ 'Luzern',
    
    grepl('Neuchatel', region) & iso.country.code == 'CH' ~ 'Neuchatel',
    
    grepl('Nidwalden', region) & iso.country.code == 'CH' ~ 'Nidwalden',
    
    grepl('Obwalden', region) & iso.country.code == 'CH' ~ 'Obwalden',
    
    grepl('Sankt Gallen', region) & iso.country.code == 'CH' ~ 'Sankt Gallen',
    
    grepl('Schaffhausen', region) & iso.country.code == 'CH' ~ 'Schaffhausen',
    
    grepl('Schwyz', region) & iso.country.code == 'CH' ~ 'Schwyz',
    
    grepl('Solothurn', region) & iso.country.code == 'CH' ~ 'Solothurn',
    
    grepl('Thurgau', region) & iso.country.code == 'CH' ~ 'Thurgau',
    
    grepl('Ticino', region) & iso.country.code == 'CH' ~ 'Ticino',
    
    grepl('Uri', region) & iso.country.code == 'CH' ~ 'Uri',
    
    grepl('Valais', region) & iso.country.code == 'CH' ~ 'Valais',
    
    grepl('Vaud', region) & iso.country.code == 'CH' ~ 'Vaud',
    
    grepl('Zug', region) & iso.country.code == 'CH' ~ 'Zug',
    
    grepl('Zurich', region) & iso.country.code == 'CH' ~ 'Zurich',
    
    
    # France - FR
    grepl('Auvergne-Rhone-Alpes', region) & iso.country.code == 'FR' ~ 'Auvergne-Rhone-Alpes',
    
    grepl('Bourgogne-Franche-Comt', region) & iso.country.code == 'FR' ~ 'Bourgogne-Franche-Comt',
    
    grepl('Bretagne', region) & iso.country.code == 'FR' ~ 'Bretagne',
    
    grepl('Centre-Val de Loire', region) & iso.country.code == 'FR' ~ 'Centre-Val de Loire',
    
    grepl('Corse', region) & iso.country.code == 'FR' ~ 'Corse',
    
    grepl('Grand-Est', region) & iso.country.code == 'FR' ~ 'Grand-Est',
    
    grepl('Hauts-de-France', region) & iso.country.code == 'FR' ~ 'Hauts-de-France',
    
    grepl('Ile-de-France', region) & iso.country.code == 'FR' ~ 'Ile-de-France',
    
    grepl('Normandie', region) & iso.country.code == 'FR' ~ 'Normandie',
    
    grepl('Nouvelle-Aquitaine', region) & iso.country.code == 'FR' ~ 'Nouvelle-Aquitaine',
    
    grepl('Occitanie', region) & iso.country.code == 'FR' ~ 'Occitanie',
    
    grepl('Pays-de-la-Loire', region) & iso.country.code == 'FR' ~ 'Pays-de-la-Loire',
    
    grepl("Provence-Alpes-Cote-d'Azur|Sainte Croix", region) & iso.country.code == 'FR' ~ "Provence-Alpes-Cote-d'Azur",

    # Mopping up
    grepl('Jinja', region) ~ 'Jinja',
    is.na(region) ~ NA,
    grepl('Devexer', region) ~ NA, 
    .default = NA
    
  )) %>% 
    left_join(iso2, by = join_by(iso.country.code, iso.subdivision.name)) %>%
  
  
  # section of genome analysed
  mutate(full.length = case_when(seqlength > 9000 ~ 'nflg',
                                 .default = 'partial')) %>%
  rename(isolation.source = isolation_source) %>%
  rename(seq.length = seqlength) %>%
  rename(seq.technology = sequencing.technology) %>%
  relocate(
    accession,
    organism,
    isolate,
    collection.date.formatted,
    collection.date.decimal,
    lab.passaged,
    iso.country.code, 
    iso.country,
    iso.subdivision.code,
    iso.subdivision.name,
    genotype,
    host.class,
    host.order,
    host.genus,
    host.species, 
    isolation.source, 
    seq.length,
    full.length,
    seq.technology,
    title,
    pubmed.id_1,
    citation_1) %>%
  mutate(genbank.createdate = parse_date(createdate)) %>%
  select(-c(createdate, updatedate)) %>% 
  mutate(is.europe = case_when(iso.country %in% c('Senegal', 'Uganda', 'South Africa', 'Central African Republic', 'Israel') ~ 'Non-European',
                               .default = 'European')) %>%
  mutate(segment = case_when(grepl('NS5|non-structural protein 5 gene', title) ~ 'NS5', 
                             grepl('envelope|E protein', title) ~ 'envelope', 
                             grepl('nflg', full.length)~'nflg')) 

  

filtered_data <- formatted_data %>%
  # only USUV 
  filter(organism == 'Usutu virus') %>%
  
  # Non passaged
  filter(lab.passaged == 0) %>%
  select(-c(strain, organism, lab.passaged)) 


# Parse NCBI Virus data


flagged_data <- filtered_data %>%
  mutate(flags_date = case_when(
    is.na(collection.date.formatted) ~ 'collection.date_missing',
    collection.date.decimal %%1==0 ~ 'collection.date_resolution',
    .default = NA)) %>%
  mutate(flags_country = case_when(
    is.na(iso.country.code) ~ 'country_missing',
    .default = NA)) %>%
  mutate(flags_subdivision = case_when(
    is.na(iso.subdivision.code) ~ 'subdivision_missing',
    .default = NA)) %>%
  mutate(flags_host = case_when(
    is.na(host.order)  ~ 'host_missing',
    .default = NA)) %>%
  unite(flags, contains('flags'), sep = ',', na.rm = T)
  

write.csv(flagged_data, './data/genbank_usuv_20231206.csv')

 


###################################################################################################