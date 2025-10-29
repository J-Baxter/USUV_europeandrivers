################################################################################
## Script Name:        Check spatial uncertainty
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-08-20
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(sf)

#Note this function assumes a lat-long ordering for coordinates
ReadKML <- function(kml_filepath, crs =  'WGS84', reverse_coords){
  kml <- scan(file = kml_filepath,
             what="",
             sep="\n",
             quiet=T)
  
  # Number of polygons present in KML
  n <- str_count(kml, 'polygon id=') %>% sum()
  poly <- c()
  
  if(n<1){
    stop('File does not contain any polygons.\n\n')
    
    }else if(n==1){
      cat('Single polygon present.')
      poly <- str_extract(kml, '-{0,1}[0-9]+\\.[0-9]+,-{0,1}[0-9]+\\.[0-9]+') %>%
      discard(is.na) %>%
      str_split(., ',') %>%
      lapply(., as.numeric) %>%
      do.call(rbind, .) %>%
      {if(isTRUE(reverse_coords)) .[,2:1] else .} %>%
      list(.) %>%
      st_polygon() %>%
      st_sfc(.) %>%
      st_set_crs(., crs)
      
  }else{
    cat('Multiple polygons present.\n\n')
    poly <- str_extract(kml, '-{0,1}[0-9]+\\.[0-9]+,-{0,1}[0-9]+\\.[0-9]+') %>%
      split(cumsum(is.na(.))) %>%     
      lapply(., function(g) g[!is.na(g)]) %>%
      compact() %>%
      lapply(., str_split, ',') %>%
      lapply(., function(x) lapply(x, as.numeric)) %>%
      lapply(., function(x) do.call(rbind, x)) %>%
      lapply(., function(x) if(isTRUE(reverse_coords)) x[,2:1] else x) %>%
      lapply(., list) %>%
      lapply(., st_polygon) %>%
      st_sfc() %>%
      st_set_crs(., crs)
  }
  
  return(poly)
}


################################### DATA #######################################
# Read and inspect data
kml_files <- list.files('./2025Jun24/kmls',
                        pattern = 'kml$',
                        full.names = T)
 
#matches <- metadata_with_concat %>%
  #filter(tipnames %in% temp_tree$tip.label) %>%
  #mutate(seq_id = coalesce(sequence_accession, sequence_isolate) %>%
           #gsub('\\/', '_', .)) %>%
  #pull(seq_id) %>%
  #paste(collapse = '|')


#kml_files <- kml_files[grepl(matches, kml_files)]


kml_poly <- lapply(kml_files, 
                   ReadKML, 
                   reverse_coords = TRUE)


metadata_with_concat <- read_csv('./data/USUV_metadata_2025Jun24_withconcatenated.csv')


################################### MAIN #######################################
# Main analysis or transformation steps
point_list <- metadata_with_concat %>% 
  
  filter(tipnames %in% temp_tree$tip.label) %>%

  # format unique IDs
  mutate(seq_id = coalesce(sequence_accession, sequence_isolate) %>%
           gsub('\\/', '_', .)) %>%
  
  # Filter seq ids present in selected KML files
  filter(seq_id %in% str_extract(kml_files,
                                 "(?<=kmls\\/).*(?=\\.kml)")) %>%
  mutate(seq_id = factor(seq_id, levels = str_extract(kml_files, 
                                                      "(?<=kmls\\/).*(?=\\.kml)"))) %>%
  arrange(seq_id) %>%
  
  # make as spatial point
  dplyr::select(seq_id, geocode_coords) %>%
  mutate(geocode_coords =  gsub('c\\(|\\)', '', geocode_coords)) %>%
  separate(geocode_coords, 
           into = c("lon", "lat"), 
           sep = ", ",
           convert = TRUE) %>%
  st_as_sf(coords = c('lon','lat'), 
           crs = 'WGS84',
           sf_column_name = 'coords') %>%
  group_split(seq_id)

# For each seq_id, is the centroid definitively within the boundaries of a 
# specified polygon?
point_within <- mapply(st_within,
                       point_list,
                       kml_poly, 
                       sparse = FALSE,
                       SIMPLIFY = F) %>%
  #unlist() %>%
  setNames(str_extract(kml_files, "(?<=kmls\\/).*(?=\\.kml)")) %>%
  enframe(name = 'seq_id',
          value = 'within_polygon') %>%
  # Check across all possible polygons
  mutate(within_polygon = lapply(within_polygon, as.vector)) %>%
  unnest_longer(within_polygon) %>%
  summarise(within_polygon = any(within_polygon),
            .by = seq_id)

# Summarise
point_within %>%
  count(within_polygon)

#point_within %>% 
  #rowid_to_column() %>%
  #filter(within_polygon == FALSE) %>% 
  #pull(rowid) %>%
  #lapply(., function(x) ggplot() + 
           #geom_sf(data = kml_poly[[x]]) + 
           #geom_sf(data = point_list[[x]]))

# Individual tests
#st_within( kml_poly[[3]],point_list[[3]], sparse = FALSE)


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################