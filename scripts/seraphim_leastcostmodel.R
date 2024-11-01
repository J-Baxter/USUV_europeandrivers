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
localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/D/"
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

features <- colnames(env_data)[c(5,8,9,10,26,30,33,34,35,36)]
"popcount"        "gcrop"           "gpast"           "gurbn"           "elevation"       "Flyway_Apod"     "livestock"      
 "Flyway_Ans"      "Flyway_Pas"      "wetland_combine"
############################################## MAIN ################################################


LeastCostModel <- function(rast, var_name, file_prefix, localTreesDirectory, nberOfExtractionFiles){
  names(rast) = var_name
  
  rast[rast[]<0] = 0
  rast[] = rast[] + 1
  
  envVariables = list(rast)
  pathModel = 2 #least-cost
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
  
  return(results)
  
}


D_env <- mapply(LeastCostModel,
       raster_list[c(5,8,9,10,26,30,33,34,35,36)-2],
       features,
       file_prefix = 'USUV_D',
       localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/D/",
       nberOfExtractionFiles= 100,
       SIMPLIFY = FALSE)

D_env %<>% bind_rows()

A_env <- mapply(LeastCostModel,
                raster_list[c(5,8,9,10,26,30,33,34,35,36)-2],
                features,
                file_prefix = 'USUV_A',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/A/",
                nberOfExtractionFiles= 100,
                SIMPLIFY = FALSE) %>%
  bind_rows()


B_env <- mapply(LeastCostModel,
                raster_list[c(5,8,9,10,26,30,33,34,35,36)-2],
                features,
                file_prefix = 'USUV_B',
                localTreesDirectory = "./2024Oct20/alignments/concatenated_alignments/B/",
                nberOfExtractionFiles= 100,
                SIMPLIFY = FALSE) %>%
  bind_rows()


combined_output <- bind_rows(A_env,
                             B_env,
                             D_env)

write_csv(combined_output, './2024Oct20/alignments/concatenated_alignments/leastcostmodelresults.csv')


#plotRaster(rast, addAxes=TRUE, addLegend=TRUE)








############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################