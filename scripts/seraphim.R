# Phylogeography visualisation
library(seraphim)
library(giscoR)
library(diagram)

# seraphim plots

allTrees <- scan(file = './2024Oct20/alignments/subset_alignments/run/USUV_2024Oct20_NFLG_A_subsampled_traits_1000.trees',
                 what = '',
                 sep = '\n',
                 quiet = T)
localTreesDirectory = "./2024Oct20/alignments/subset_alignments/A/"
burnIn <- 0
randomSampling <- FALSE
nberOfTreesToSample <- 1000
mostRecentSamplingDatum <- decimal_date(c_most_recent_date %>% ym())
coordinateAttributeName <- "location"
treeExtractions(localTreesDirectory,
                allTrees,
                burnIn, 
                randomSampling, 
                nberOfTreesToSample, 
                mostRecentSamplingDatum,
                coordinateAttributeName)


mcc_tree <-readAnnotatedNexus('./2024Oct20/alignments/subset_alignments/run/USUV_2024Oct20_NFLG_A_subsampled_traits_mcc.tree')
mcc_tree_tbl <- as_tibble(mcc_tree)
source("./imported_scripts/mccExtractions.r") # Script obtained from the GitHub tutorial folder.

mcc_tab <- mccExtractions(mcc_tre, mostRecentSamplingDatum)

# Step 4: Estimating the HPD region for each time slice ----
nberOfExtractionFiles <- nberOfTreesToSample
prob <- 0.95
precision <- 0.1 # time interval that will be used to define the successive time slices
startDatum <- min(mcc_tab[,"startYear"])
polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

# Step 6: Co-plotting the HPD regions and MCC tree ----
# Step 6: Co-plotting the HPD regions and MCC tree ----
map <-all_europe_sf
# plot nodes on map
nodes <- mcc_tre_tbl %>%
  dplyr::select(node,height, location1, location2) %>%
  mutate(year = decimal_date(ym(c_most_recent_date)) - as.numeric(height)) %>%
  st_as_sf(coords = c( 'location2', 'location1'), 
           crs = 4326)


arrows <-  mcc_tre_tbl %>%
  dplyr::select(node, parent) %>%
  left_join(mcc_tre_tbl %>%
              dplyr::select(node, location1, location2)) %>%
  left_join(mcc_tre_tbl %>%
              dplyr::select(node, location1, location2),
            by = join_by(parent== node)) %>%
  mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
  mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                 .default = location1.y),
         location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                 .default = location2.y))

sf_list <- lapply(polygons, st_as_sf) 
sf_combined <- bind_rows(sf_list) %>%
  pivot_longer(cols = starts_with('2'), names_to = 'year', values_to = 'value') %>%
  filter(value == 1) %>%
  mutate(year = as.numeric(year)) %>%
  dplyr::select(-value) %>%
  st_set_crs(4326) 


ggplot(map) +
  geom_sf() +
  
  # Plot HPD polygons
  geom_sf(data = sf_combined ,aes(fill = year),lwd = 0,
          alpha = 0.2) + 
  
  
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
  
  
  coord_sf(ylim = c(40,60), xlim = c(-5, 20), expand = FALSE) +
  theme_void(base_size = 18) + 
  guides(fill = guide_colourbar(
                                theme = theme(
                                  legend.key.width  = unit(15, "lines"),
                                  legend.key.height = unit(1, "lines")
                                  ))) + 
  
  theme(legend.background = element_rect(colour="white", fill="white"),
        legend.position =c(.7,.1),
        legend.title = element_text( vjust = 1),
        legend.direction="horizontal")




# Phylogeography dispersal vectors

nberOfExtractionFiles = 100
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = FALSE
outputName = "USUV_C"
nberOfCores = 1
slidingWindow = 1

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, 
                 showingPlots, outputName, nberOfCores, slidingWindow)