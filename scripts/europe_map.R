# paper figures
library(terra)
library(sf)
library(rnaturalearth)
library(giscoR)
library(eurostat)
library(tidyverse)


# Initialise map of europe subdivisions (NUTS level 2)
# Data obtained from Eurostat GISCO (the Geographic Information System of the COmmission) via API 
# package giscoR

nuts2 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326", #WGS84 projection
  resolution = "10", #1:10million
  nuts_level = "2"
)

ggplot(nuts2) + 
  geom_sf() + 
  coord_sf(ylim = c(34,72), xlim = c(-12, 45), expand = FALSE)
