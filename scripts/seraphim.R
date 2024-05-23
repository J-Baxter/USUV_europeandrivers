# Phylogeography visualisation
library(seraphim)

# seraphim plots

allTrees <- scan(file = '~/Downloads/USUV_EU3_nflg_2024May7_subsampled.trees',
                 what = '',
                 sep = '\n',
                 quiet = T)
localTreesDirectory = "./Extracted_trees/"
burnIn <- 0
randomSampling <- FALSE
nberOfTreesToSample <- 1000
mostRecentSamplingDatum <- 2021
coordinateAttributeName <- "location"
treeExtractions(localTreesDirectory,
                allTrees,
                burnIn, 
                randomSampling, 
                nberOfTreesToSample, 
                mostRecentSamplingDatum,
                coordinateAttributeName)


mcc_tre <-readAnnotatedNexus('~/Downloads/USUV_EU3_nflg_2024May7_subsampled_mcc.tree')
source("~/Downloads/mccExtractions.R") # Script obtained from the GitHub tutorial folder.

mcc_tab <- mccExtractions(mcc_tre, mostRecentSamplingDatum)

# Step 4: Estimating the HPD region for each time slice ----
nberOfExtractionFiles <- nberOfTreesToSample
prob <- 0.95
precision <- 0.1 # time interval that will be used to define the successive time slices
startDatum <- min(mcc_tab[,"startYear"])
polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


# Step 5: Defining the different colour scales to use ----
colour_scale <- colorRampPalette(brewer.pal(9,"YlGnBu"))(141)[21:121]
minYear <- min(mcc_tab[,"startYear"])
maxYear <- max(mcc_tab[,"endYear"])
endYears_indices <- (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours <- colour_scale[endYears_indices]
polygons_colours <- rep(NA, length(polygons))
for (i in 1:length(polygons)){
  date <- as.numeric(names(polygons[[i]]))
  polygon_index <- round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] <- paste0(colour_scale[polygon_index],"40")
}

# Step 6: Co-plotting the HPD regions and MCC tree ----
# Step 6: Co-plotting the HPD regions and MCC tree ----
map <- ne_countries(returnclass = "sf") %>% st_transform(., 4326) filter(COUNTRY%in% c(
  'United Kingdom','France','Spain', 'Portugal', 'Andorra', 'Ireland', 'Poland',
  'Monaco', 'Switzerland', 'Italy', 'Luxembourg', 'Netherlands', 'Moldova', "Bosnia and Herzegovina",
  'Belgium', 'Netherlands', 'Germany', 'Denmark', 'Norway', 
  'Sweden', 'Finland', 'Estonia', 'Lithuania', 'Latvia', 
  'Belarus', 'Ukraine', 'Turkey', 'Greece', 'Cyprus', "North Macedonia", "Liechtenstein",
  'Slovakia', 'Slovenia', 'Croatia', 'Serbia', 'Hungary', 'Romania', 'Bulgaria', 
  'Czechia', 'Austria', 'Albania', 'Kosovo', 'Guernsey', 'Jersey', 'Isle of Man', 'Faroe Islands',
  'Iceland', 'Monaco', 'San Marino', 'Russia', 'Morocco', 'Algeria', 'Tunisia', "Libya", 'Egypt', 'Israel', 'Syria', 'Lebanon'))%>% 
   # Set CRS to lat-long.


#%>%
#st_as_sf(coords = c("location1", "location2"), crs = 4326) %>% 
 # st_transform(crs = crs_use)
# Plot the continuous phylogeography in an external PDF:
pdf("test.pdf", width = 6, height = 6.3)
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")

plot(st_geometry(map),
     col="grey", border = "#D1D1D1", lwd = 1.5,
     xlim = c(4, 25),
     ylim = c(45, 65))
for (i in 1:length(polygons)){
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tab)[1]){
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc_tab)[1]:1){
  if (i == 1)
  {
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16, col=colour_scale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1, col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8)
}

# 2nd add the axes:
#axis(side = 1, 
    # at = seq(-5, -3.8, 0.1), 
    # pos= 55.372, # originally at 55.4
    # mgp=c(0,0.2,0), 
     #cex.axis=0.5, 
     #lwd=1, 
     #lwd.tick=0.2, 
     #padj=-0.8, 
     #tck=-0.01, 
     #col.axis="gray30")
#axis(side = 2, 
    # at = seq(55.3, 56.1, 0.1), 
    # pos= -4.877, # originally at -5
    # mgp=c(0,0.2,0), 
     #cex.axis=0.5, 
     #lwd=1, 
     #lwd.tick=0.2, 
     #padj=-0.8, 
     #tck=-0.01, 
     #col.axis="gray30")

# 3rd add the legend:
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab[,"startYear"]); rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.95,0.08,0.1),
     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0), at=seq(2008,2023,1)))
dev.off()


# Phylogeography dispersal vectors

nberOfExtractionFiles = 100
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = FALSE
outputName = "WNV"
nberOfCores = 1
slidingWindow = 1

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, 
                 showingPlots, outputName, nberOfCores, slidingWindow)