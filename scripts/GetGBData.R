# Initial Data Collection
library(rentrez)
library(genbankr)
source('./scripts/apikey.R')


###################################################################################################
# Function to extract Genbank collection dates from NCBI nucleotide database
# Requires function NCBI API key, otherwise downloads will exceed maximum number of calls


# Requires pubication ID

GetGBData <- function(searchterm, key = API_KEY, max = 10){
  
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
    
  }else{
    test_summ <- split(test_id,
                       ceiling(seq_along(test_id) / 200)) %>%
      lapply(.,  entrez_summary, db="nuccore") %>%
      list_flatten() %>%
      lapply(., function(x) keep(x, is.character)) %>%
      lapply(., as_tibble) %>%
      bind_rows()%>%
      rename(accession = caption) 
    
  }
  
  
  #GENBANKR search
  accn <- test_summ %>%
    select(accession) %>%
    unlist() %>%
    as.vector()
  
  MyFunc1 <- function(accn, key = key){
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
  
  list_data <- lapply(accn, MyFunc1) %>% 
    bind_rows()
  
  combined_dataframe <- full_join(test_summ, 
                                  list_data, 
                                  by = join_by(accession))
  
  
  
  out <- combined_dataframe %>% 
    mutate(result = map(accession, getAnnotationsGenBank)) %>%  
    unnest(result) %>%  
    select(-c(type, product, others))

  return(out)
}


full_data <- GetGBData("usutu[All Fields] AND viruses[filter]", 
                  max = 2000)
full_data <- test


# Country codes
# Region detail 
# Host organism 
# Host type
# format year
# decimal year
# section of genome analysed
# sequence technology

data_formatted <- full_data %>%
  select(-c(term, flags, biomol, topology, segsetsize, sourcedb, tech, completeness, projectid, genome, contains('assembly')))





###################################################################################################