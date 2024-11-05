####################################################################################################
####################################################################################################
## Script name: Draft Seraphim Least Cost Model
##
## Purpose of script:
##
## Date created: 2024-10-31
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 4) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(seraphim)

# User functions


############################################## DATA

localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/C/"
burnIn <- 0
coordinateAttributeName <- "location"
treeExtractions(localTreesDirectory,
randomSampling <- TRUE
nberOfTreesToSample <- 500
mostRecentSamplingDatum <- decimal_date(most_recent_date%>% ymd())
                allTrees,
                burnIn, 
                randomSampling, 
                nberOfTreesToSample, 
                mostRecentSamplingDatum,
                coordinateAttributeName,
                nberOfCores = 6)


localTreesDirectory = "Extracted_trees"
allTrees = scan(file="RABV_gamma.trees", what="", sep="\n", quiet=TRUE)


burnIn <- 0
randomSampling <- TRUE
nberOfTreesToSample <- 500
mostRecentSamplingDatum <- decimal_date( %>% ymd())
coordinateAttributeName <- "location"

treeExtractions(localTreesDirectory,
                allTrees,
                burnIn, 
                randomSampling, 
                nberOfTreesToSample, 
                mostRecentSamplingDatum,
                coordinateAttributeName,
                nberOfCores = 16)


# Prepare environmental rasters
env_data  <- readRDS("./data/Predictors_tep_20231014.Rdata") %>% 
  dplyr::mutate_all(as.numeric)

rast <- env_data %>%
  select(lon, lat, wetland_combine) %>%
  rasterFromXYZ()


raster_list <- lapply(3:ncol(env_data),
                      function(i) env_data %>%
                        select(lon, lat, i) %>%
                        rasterFromXYZ())

features <- colnames(env_data)[3:ncol(env_data)]

############################################## MAIN ################################################


LeastCostModel <- function(rast, var_name, file_prefix, localTreesDirectory, nberOfExtractionFiles, pathModel){
  names(rast) = var_name
  
  rast[rast[]<0] = 0
  rast[] = rast[] + 1
  
  envVariables = list(rast)
  
  # path model values:
  # > 0 -> dispersal direction (least-cost)
  # > 2 -> dispersal velocity (least-cost)
  
  if(pathModel == 2){
    resistances = list(TRUE)
    avgResistances = list(TRUE)
    fourCells = FALSE
    nberOfRandomisations = 0
    randomProcedure = 3
    outputName = paste(file_prefix, var_name,'least-cost', sep = '_')
    
    spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables,
                  pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
                  nberOfCores = 6,
                  randomProcedure, outputName)
    
    temp <- paste(file_prefix, var_name,'least-cost_linear_regression_results.txt', sep = '_')
    print(temp)
    table = read.table(temp,
                       header=T)
    
    col <- paste0("Univariate_LR_coefficients_", var_name[1], "_R")
    
    LR_coefficients = table[,col]
    print(sum(LR_coefficients > 0))
    
    col <- paste0("Univariate_LR_delta_R2_", var_name[1], "_R")
    
    Qs = table[,col]
    print(sum(Qs > 0))
    
    nberOfRandomisations = 1
    
    spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables,
                  pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
                  randomProcedure, outputName)
    
    temp <- paste(file_prefix, var_name,'least-cost_randomisation_Bayes_factors.txt', sep = '_')
    bf <- read.table(temp, header = T)
    
    results <- tibble(var = var_name,
                      prefix = file_prefix, 
                      LR_sum = sum(LR_coefficients > 0),
                      Qs_sum = sum(Qs > 0),
                      bayes_factor = bf[2,1])
  }else if(pathModel == 0){
    
    resistances = list(FALSE)
    avgResistances = list(FALSE)
    fourCells = FALSE
    nberOfRandomisations = 1
    randomProcedure = 3
    outputName = paste(file_prefix, var_name,'least-cost', sep = '_')
    
    spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables,
                  pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
                  nberOfCores = 6,
                  randomProcedure, outputName) 
    
    temp_E <- paste(file_prefix, var_name,'least-cost_direction_E_Bayes_factors.txt', sep = '_')
    temp_R <- paste(file_prefix, var_name,'least-cost_direction_R_Bayes_factors.txt', sep = '_')
    
    bf_E <- read.table(paste0('./', temp_E), header = T) 
    bf_R <- read.table(paste0('./', temp_R), header = T) 
    
    results <- bind_rows(bf_E, bf_R) %>%
      as_tibble(., rownames = 'var')
  }
 
  
  return(results)
  
}

A_velocity <- mapply(LeastCostModel,
                raster_list[-c(1, 2, 36, 37)],
                features[-c(1, 2, 36, 37)],
                file_prefix = 'USUV_A',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/A/",
                nberOfExtractionFiles= 100,
                pathModel = 2,
                SIMPLIFY = FALSE) %>%
  bind_rows()



B_velocity <- mapply(LeastCostModel,
                     raster_list[-c(1, 2, 36, 37)],
                     features[-c(1, 2, 36, 37)],
                file_prefix = 'USUV_B',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/B/",
                nberOfExtractionFiles= 100,
                pathModel = 2,
                SIMPLIFY = FALSE) %>%
  bind_rows()

C_velocity <- mapply(LeastCostModel,
                     raster_list[-c(1, 2, 36, 37)],
                     features[-c(1, 2, 36, 37)],
                file_prefix = 'USUV_C',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/C/",
                nberOfExtractionFiles= 100,
                pathModel = 2,
                SIMPLIFY = FALSE)

C_velocity %<>% bind_rows()

D_velocity <- mapply(LeastCostModel,
                     raster_list[-c(1, 2, 36, 37)],
                     features[-c(1, 2, 36, 37)],
                file_prefix = 'USUV_D',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/D/",
                nberOfExtractionFiles= 100,
                pathModel = 2,
                SIMPLIFY = FALSE)

D_velocity %<>% bind_rows()

combined_output <- bind_rows(A_velocity,
                             B_velocity,
                             C_velocity,
                             D_velocity)

write_csv(combined_output, './2024Oct20/alignments/concatenated_alignments/leastcostmodelresults_velocity.csv')


#plotRaster(rast, addAxes=TRUE, addLegend=TRUE)

A_direction <- mapply(LeastCostModel,
                      raster_list[-c(1, 2, 36, 37)],
                      features[-c(1, 2, 36, 37)],
                      file_prefix = 'USUV_A',
                      localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/A/",
                      nberOfExtractionFiles= 100,
                      pathModel = 0,
                      SIMPLIFY = FALSE)

B_direction <- mapply(LeastCostModel,
                      raster_list[-c(1, 2, 36, 37)],
                      features[-c(1, 2, 36, 37)],
                      file_prefix = 'USUV_B',
                      localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/B/",
                      nberOfExtractionFiles= 100,
                      pathModel = 0,
                      SIMPLIFY = FALSE)

C_direction <- mapply(LeastCostModel,
                      raster_list[-c(1, 2, 36, 37)],
                      features[-c(1, 2, 36, 37)],
                      file_prefix = 'USUV_C',
                      localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/C/",
                      nberOfExtractionFiles= 100,
                      pathModel = 0,
                      SIMPLIFY = FALSE)

A_direction %<>%
  bind_rows() %>%
  drop_na(BF)
B_direction %<>%
  bind_rows() %>%
  drop_na(BF)
C_direction %<>%
  bind_rows() %>%
  drop_na(BF)

combined_output <- bind_rows(A_direction,
                             B_direction,
                             C_direction,
                             D_direction)

write_csv(combined_output, './2024Oct20/alignments/concatenated_alignments/leastcostmodelresults_direction.csv')


D_direction <- mapply(LeastCostModel,
                     raster_list[-c(1, 2, 36, 37)],
                     features[-c(1, 2, 36, 37)],
                     file_prefix = 'USUV_D',
                     localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/D/",
                     nberOfExtractionFiles= 100,
                     pathModel = 0,
                     SIMPLIFY = FALSE)

D_direction %<>%
  bind_rows() %>%
  drop_na(BF)
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################