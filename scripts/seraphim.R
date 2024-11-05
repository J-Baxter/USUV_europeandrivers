# Phylogeography visualisation
library(seraphim)
library(giscoR)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(diagram)
library(lubridate)
library(tidyverse)
library(treeio)
# seraphim plots

allTrees <- scan(file = './2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_noFLI_concat_C_subsampled_temp_1000.trees',
                 what = '',
                 sep = '\n',
                 quiet = T)
TREEFILE <- './2024Oct20/alignments/concatenated_alignments/USUV_2024Oct20_noFLI_concat_C_subsampled_temp_mcc.tree'
mcc_tree <-read.beast(TREEFILE)
mcc_tree_tbl <- as_tibble(mcc_tree)



metadata_in_tree <- read_csv('./data/USUV_metadata_noFLI_2024Oct20_withconcatenated.csv') %>%
  dplyr::filter(tipnames %in% mcc_tree@phylo$tip.label) 

stopifnot(nrow(metadata_in_tree) == Ntip(mcc_tree@phylo)) #sanity check

most_recent_date <- metadata_in_tree %>%
  pull(date_ymd) %>% #note explicit assumption that most recent date will be ymd not ym
  max(na.rm = TRUE)


localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/C/"
burnIn <- 0
randomSampling <- TRUE
nberOfTreesToSample <- 500
mostRecentSamplingDatum <- decimal_date(most_recent_date%>% ymd())
coordinateAttributeName <- "location"
treeExtractions(localTreesDirectory,
                allTrees,
                burnIn, 
                randomSampling, 
                nberOfTreesToSample, 
                mostRecentSamplingDatum,
                coordinateAttributeName,
                nberOfCores = 6)


source('./scripts/import_scripts_2/mccExtractions.r') # Script obtained from the GitHub tutorial folder.

mcc_tab <- mccExtractions(readAnnotatedNexus(TREEFILE), mostRecentSamplingDatum)

# Step 4: Estimating the HPD region for each time slice ----
nberOfExtractionFiles <- nberOfTreesToSample
prob <- 0.95
precision <- 0.1 # time interval that will be used to define the successive time slices
startDatum <- min(mcc_tab[,"startYear"])

# Format Polygons
polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, 
                                            nberOfExtractionFiles,
                                            prob, 
                                            startDatum, 
                                            precision))

sf_list <- lapply(polygons, st_as_sf) 

sf_combined <- bind_rows(sf_list) %>%
  pivot_longer(cols = starts_with('2'), names_to = 'year', values_to = 'value') %>%
  filter(value == 1) %>%
  mutate(year = as.numeric(year)) %>%
  dplyr::select(-value) %>%
  st_set_crs(4326) 

# Set Base Map
nongisco_shapefile <- ne_countries(scale = 10, country = c('bosnia and herzegovina', 'kosovo', 'andorra'), returnclass = "sf") %>%
  dplyr::select(name_en,
                iso_a2_eh,
                geometry) %>%
  rename(CNTR_CODE = iso_a2_eh,
         NAME_LATN = name_en) %>%
  mutate(NUTS_ID = CNTR_CODE)

nuts1 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "03", #1:10million
  nuts_level = "1")

all_europe_sf <- bind_rows(nuts1, nongisco_shapefile)

map <-all_europe_sf


# Format Nodes
nodes <- mcc_tree_tbl %>%
  dplyr::select(node,height, location1, location2) %>%
  mutate(year = decimal_date(ymd(most_recent_date)) - as.numeric(height)) %>%
  st_as_sf(coords = c( 'location2', 'location1'), 
           crs = 4326)

# Format Arrows
arrows <-  mcc_tree_tbl %>%
  dplyr::select(node, parent) %>%
  left_join(mcc_tree_tbl %>%
              dplyr::select(node, location1, location2)) %>%
  left_join(mcc_tree_tbl %>%
              dplyr::select(node, location1, location2),
            by = join_by(parent== node)) %>%
  mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
  mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                 .default = location1.y),
         location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                 .default = location2.y))

# Plot in GGplot
ggplot(map) +
  geom_sf() +
  
  # Plot HPD polygons
  geom_sf(data = sf_combined ,aes(fill = year),lwd = 0,
          alpha = 0.05) + 
  
  
  # Plot Branches (support level?)
  geom_curve(data = arrows, aes(x = location2.x, y = location1.x,
                                xend = location2.y, yend = location1.y),
             lwd = 0.3,
             curvature = 0.2) + 
  
  # Plot nodes
  geom_sf(data = nodes, shape = 21 , size = 2, aes(fill = year))+
  
  # Set colour scale
  scale_fill_viridis_c('Year',
                       direction = -1,
                       option = 'C')+
  
  
  coord_sf(ylim = c(35,60), xlim = c(-8, 33), expand = FALSE) +
  theme_void(base_size = 18) + 
  guides(fill = guide_colourbar( )) + 
  
  theme(legend.background = element_rect(colour="white", fill="white"),
        legend.position =c(.7,.1),
        legend.title = element_text( vjust = 1),
        legend.direction="horizontal",
        legend.key.width  = unit(3, "lines"),
        legend.key.height = unit(1, "lines"))




# Phylogeography dispersal vectors
localTreesDirectory <-"./2024Oct20/alignments/concatenated_alignments/C/"
nberOfExtractionFiles = 100
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = FALSE
outputName = "USUV_C"
nberOfCores = 6
slidingWindow = 1

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, 
                 showingPlots, outputName, nberOfCores, slidingWindow)