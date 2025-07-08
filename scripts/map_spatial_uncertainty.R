################################################################################
## Script Name:        Tip Location Sampling Uncertainty
## Purpose:            Edit Beauti XMLs to incorporate location uncertainty 
## Author:             James Baxter
## Date Created:       2025-07-07
## To-do: 1) VectorNet with aggregated mosquito probabilities
##        2) connection between KML filenames and sequences - /
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
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(giscoR)

# Function adpated from https://github.com/sdellicour/h5n1_mekong
WriteNUTS3KML <- function(x, prefix = NULL){
  seq_id <- x$tipnames
  sampling_probability = x$sampling_probability
  coords <- x$coords[[1]]
  
  sink(file=paste(prefix, seq_id, ".kml", sep=""))
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  cat("\n")
  cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">")
  cat("\n")
  cat(paste("\t<polygon id=\"", paste(seq_id,sep="_"), 
            "\" samplingProbability=\"", 1, "\">", sep=""))
  cat("\n")
  cat("\t\t<coordinates>")
  cat("\n")
  
  for (j in 1:dim(coords)[1]) {
    cat(paste("\t\t\t", coords[j, 2], ",", coords[j, 1], ",0", sep=""))
    cat("\n")
  }
  cat("\t\t</coordinates>")
  cat("\n")
  cat("\t</polygon>")
  cat("\n")
  cat("</kml>")
  cat("\n")
  sink(NULL)
  closeAllConnections()
}


# Function adpated from https://github.com/sdellicour/h5n1_mekong
WriteHeteroKML <- function(x, prefix = NULL){

}


# Function from https://github.com/sdellicour/h5n1_mekong
UpdateBEAUti <- function(xml_filepath){
  xml = scan(file="BEAST_template.xml", what="", sep="\n", quiet=T)
  directory = "H5N1_polygons"
  sink(file="H5N1_clade1.xml")
  for (i in 1:length(xml)) {
    cat(xml[i]); cat("\n")
    if (xml[i]=="\t</continuousDiffusionStatistic>") {
      cat("\n")
      for (j in 1:length(sequenceIDs)) {
        cat(paste("\t<leafTraitParameter id=\"",sequenceIDs[j],".trait\" taxon=\"",names[j],"\">",sep="")); cat("\n")
        cat(paste("\t\t<treeModel idref=\"treeModel\"/>",sep="")); cat("\n")
        cat(paste("\t\t<parameter idref=\"leaf.location\"/>",sep="")); cat("\n")
        cat(paste("\t</leafTraitParameter>",sep="")); cat("\n")
      }
      cat("\n")
      for (j in 1:length(sequenceIDs)) {
        cat(paste("\t<flatGeoSpatialPrior id=\"",sequenceIDs[j],"_polygons\" taxon=\"",names[j],"\" kmlFileName=\"",directory,"/",sequenceIDs[j],".kml\" inside=\"true\" union=\"true\" cache=\"true\">",sep="")); cat("\n")
        cat(paste("\t\t<data>",sep="")); cat("\n")
        cat(paste("\t\t\t<parameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
        cat(paste("\t\t</data>",sep="")); cat("\n")
        cat(paste("\t</flatGeoSpatialPrior>",sep="")); cat("\n")
      }
      cat("\n")		
    }
    if (xml[i]=="\t\t</precisionGibbsOperator>") {
      cat("\n")
      for (j in 1:length(sequenceIDs)) {
        cat(paste("\t\t<uniformGeoSpatialOperator weight=\"0.01\">",sep="")); cat("\n")
        cat(paste("\t\t\t<parameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
        cat(paste("\t\t\t<flatGeoSpatialPrior idref=\"",sequenceIDs[j],"_polygons\"/>",sep="")); cat("\n")
        cat(paste("\t\t</uniformGeoSpatialOperator>",sep="")); cat("\n")
      }
      cat("\n")
    }
    if (xml[i]=="\t\t\t\t<multivariateWishartPrior idref=\"location.precisionPrior\"/>") {
      cat("\n")
      cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
      for (j in 1:length(sequenceIDs)) {
        cat(paste("\t\t\t\t<flatGeoSpatialPrior idref=\"",sequenceIDs[j],"_polygons\"/>",sep="")) 
        cat("\n")
      }
      cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
      cat("\n")
    }					
    if (xml[i]=="\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood\"/>") {
      if (xml[i-3]=="\t\t\t<strictClockBranchRates idref=\"branchRates\"/>") {
        cat("\n")
        for (j in 1:length(sequenceIDs)) {
          cat(paste("\t\t\t<leafTraitParameter idref=\"",sequenceIDs[j],".trait\"/>",sep="")); cat("\n")
        }
        cat("\n")
      }
    }	
  }
  sink(NULL)
}


################################### DATA #######################################
# Read and inspect data
vector_net <- read_sf('./spatial_data/VectornetMAPforMOODjan21.shp', crs = st_crs(nuts0)) %>%
  st_make_valid() %>%
  st_transform(st_crs(vector_net))

metadata <- read_csv('./data/USUV_metadata_all_2025Jun24.csv')

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

nuts_all <- bind_rows(nuts0,
                      nuts1,
                      nuts2)

culex_abundance_files <- list.files('./culex_models/statistical_models/Abundace_model_predictions',
                                    full.names = T)

################################### MAIN #######################################
# Main analysis or transformation steps
# Pre-format vectornet polygons
vect_id <- as_tibble(vector_net[-1083,]) %>%
  rowid_to_column() %>%
  st_as_sf()

# Prepare predictions from abundance data
culex_abundance <- lapply(culex_abundance_files, read_csv) %>%
  setNames(gsub('\\.\\/culex_models\\/statistical_models\\/Abundace_model_predictions\\/culex_|\\.csv', '', culex_abundance_files)) %>%
  bind_rows(., .id = 'mm_yyyy') %>%
  separate_wider_delim(., mm_yyyy, '_', names = c('month', 'year')) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) 

# Simple average (mean) over all month-year combinations, aggregated by 
# polygons
culex_abundance_by_vectornet <- aggregate(culex_abundance %>% 
                                         dplyr::select(prediction,geometry),
                                         vect_id, 
                                         mean) %>%
  rowid_to_column() 

#ggplot(culex_abundance_aggregate,
      # aes(fill = prediction)) + 
 # geom_sf() + 
 # coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  #scale_fill_distiller(palette = 'RdYlGn', 
                       #transform = 'log10', 
                       #direction = -1, , 
                       #na.value="lightgrey") +
  #facet_grid(rows = vars(year), cols = vars(month)) + 
  #theme_void() 


##### NUTS3 (homogeneous sampling within NUTS3) #####
metadata %>% 
  filter(is_europe == '1') %>%
  filter(location_precision == 'nuts3') %>%
  dplyr::select(tipnames, eurostat_polygon, 'nuts3_name') %>%
  
  # Join polygons and convert to matrix form. Selecting best resolution by default
  left_join(vect_id, by = join_by('eurostat_polygon' == 'rowid')) %>%
  mutate(coords = list(unlist(geometry, recursive = FALSE)[[1]][[1]]),
         sampling_probability = 1) %>%
  dplyr::select(1,5,6) %>%
  
  # splt dataframe by sequence
  rowid_to_column(var = 'id') %>%
  group_split(id) %>%
  
  # write KML file for each sequence
  lapply(., WriteNUTS3KML, prefix = './2025Jun24/kmls/')


##### NUTS0, NUTS1, and NUTS2 (sampling determined by culex abundance ) #####
testnuts1 <- metadata %>% 
  filter(is_europe == '1') %>%
  filter(location_precision %in% c('nuts0', 'nuts1', 'nuts2')) %>%
  .[1:5,] %>%
  dplyr::select(tipnames, ends_with('_id')) %>%
  pivot_longer(cols = -1, names_to = 'nuts_level', values_to = 'nuts_id') %>%
  drop_na() %>%
  mutate(nuts_level = gsub('_id', '', nuts_level)) %>%
  
  # Join current polygon
  left_join(nuts_all %>% dplyr::select(NUTS_ID, geometry),
            by = join_by('nuts_id' == 'NUTS_ID')) %>%
  
  # list all vectornet polygons within current polygon
  mutate(eurostat_polygon = map(geometry, ~ {
    which(st_overlaps(vect_id$geometry, .x, sparse = FALSE)[, 1] | st_within(vect_id$geometry, .x, sparse = FALSE)[, 1] )
  })) %>%
  dplyr::select(-geometry) %>%
  
  # Expand dataframe so that there is 1 row per eurostate polygon
  unnest(eurostat_polygon) %>%
  
  # join eurostat polygons, with aggregated culex abundance. Selecting best resolution by default
  left_join(culex_abundance_by_vectornet,
            by = join_by('eurostat_polygon' == 'rowid')) %>%
  mutate(coords = list(unlist(geometry, recursive = FALSE)[[1]][[1]])) %>%
  
  # Calculate sampling probabilities 
  group_by(tipnames) %>%
  mutate(sampling_probability = prediction / sum(prediction, na.rm = TRUE)) %>%
  ungroup() %>%
  
  dplyr::select(1,7, 8) %>% 
  
  # splt dataframe by sequence
  group_split(tipnames) %>%
  
  # write KML file for each sequence
  lapply(., WriteHeteroKML, prefix = './2025Jun24/kmls/')



map_metadata <- metadata %>%
  
  # format geo
  mutate(geocode_coords =  gsub('c\\(|\\)', '', geocode_coords)) %>%
  separate(geocode_coords, into = c("lat", "lon"), sep = ", ", convert = TRUE) %>%
  st_as_sf(coords = c('lat','lon'), crs = st_crs(vector_net), sf_column_name = 'coords') %>%
  st_join(.,
          vect_id, 
          join = st_within, 
          # left = FALSE,
          # largest = TRUE
  ) %>% 
  drop_na(nuts2_id) %>%
  
  # format 
  #left_join(lineage_tbl, by = join_by('tipnames' == 'label')) %>%
  #drop_na(`GRI Lineage Level 1`) %>%
  count( rowid) %>%
  #rename(lineage = `GRI Lineage Level 1`) %>%
  st_drop_geometry()


map_data <- expand_grid(rowid = vect_id$rowid,
                        n = NA_integer_) %>%
  rows_patch(map_metadata,
             unmatched = 'ignore',
             by = c('rowid')) %>%
  left_join(vect_id) #%>%
  #filter(grepl('^A', lineage)) 
# filter(CNTR_CODE != 'TR') 


map_data %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n), colour = 'white') + 
  #geom_sf(colour = 'white', data = nuts0, linewidth = 0.5, alpha = 0.01) + 
  coord_sf(ylim = c(34,72), xlim = c(-11, 34), expand = FALSE) +
  scale_fill_distiller(palette = 'RdPu', 
                       #transform = 'log10', 
                       direction = -1,  
                       na.value="lightgrey") +
  #facet_wrap(~lineage) +
  theme_void() + 
  theme(legend.position = 'none')


# Exact = Exact






WriteNUTS3KML(test_nuts3[[1]])





################################### OUTPUT #####################################

# Update Beauti 




#################################### END #######################################
################################################################################