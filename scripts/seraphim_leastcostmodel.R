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


############################################## DATA ################################################
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




############################################## MAIN ################################################


rast = raster("Elevation_raster.asc")
names(rast) = "elevation"
rast[rast[]<0] = 0
rast[] = rast[] + 1

plotRaster(rast, addAxes=TRUE, addLegend=TRUE)


envVariables = list(rast)
pathModel = 2
resistances = list(TRUE)
avgResistances = list(TRUE)
fourCells = FALSE
nberOfRandomisations = 0
randomProcedure = 3
outputName = "RABV_elevation_least-cost"


spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables,
              4
              pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
              randomProcedure, outputName)

tab = read.table("RABV_elevation_least-cost_linear_regression_results.txt",
                 header=T)
LR_coefficients = tab[,"Univariate_LR_coefficients_elevation_R"]
print(sum(LR_coefficients > 0))

Qs = tab[,"Univariate_LR_delta_R2_elevation_R"]
print(sum(Qs > 0))

spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables,
              pathModel, resistances, avgResistances, fourCells, nberOfRandomisations,
              randomProcedure, outputName)
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################