################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-10-30
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


# Read in logfile and sample mean evolutionary rates
ResampleEvoRate <- function(log_filepath, 
                            burnin = 0.1, 
                            n = 100){
  
  out <- read_delim(log_filepath, 
                    delim ='\t', 
                    skip = 4) %>%
    
    # remove burnin
    filter(state >= 0.1*max(state)) %>%
    
    #select mean rate
    dplyr::select(meanRate) %>%
    
    # down sample at random - confirm this is valid. Downsample or thin?
    sample_n(n)
  
  return(out)
}


# Read in XMLs and update with fixed rate from sampled rates extracted in 
# previous step.
AdjustFixedRate <- function(xml_filepath,
                            new_rate){
  require(scales)
  
  current_xml <- xml_filepath
  
  # name iteration (i.e check if xml with rateX.xml is present and add +)
  working_directory <- gsub('[\\/](?=[^\\/]*$).*', '', current_xml, perl = T)
  xml_filestem <- gsub('.*[\\/](?=[^\\/]*$)', '', current_xml, perl = T) %>%
    gsub('.xml$', '', .)
  
  existing_xmls <- list.files(path = working_directory,
                              pattern = xml_filestem) %>%
    grepl('rate.*.xml$', .)
  
  
  if(any(existing_xmls)){
    n <- sum(existing_xmls, na.rm = TRUE) %>%
      add(1) %>%
      str_pad(., pad = '0', 3)
    
    updated_xml <- gsub('.xml$', paste0('_rate', n, '.xml') , current_xml)
    
  }else{
    updated_xml <- gsub('.xml$', '_rate001.xml', current_xml)
  }
  
  
  # Read Existing XML and extract sequence names
  xml <-  scan(file = current_xml, what="", sep="\n", quiet=T)
  
  new_rate_formatted <- scientific(new_rate, digits =8)
  
  # Write Updated XML
  sink(file = updated_xml)
  
  for (i in 1:length(xml)) {
    if(grepl('parameter id="clock.rate', xml[i])){
      xml[i] <- paste("\t\t\t<parameter id=\"clock.rate\" value=\"", new_rate_formatted,"\"/>",sep="")
    }
    cat(xml[i]); cat("\n")
  }
  
  sink(NULL)
  
}


################################### DATA #######################################
# Read and inspect data
beast_logs <- 
updated_xmls <- 


################################### MAIN #######################################
# Main analysis or transformation steps

ResampleEvoRate('./2025Jun24/europe_clusters/NFLG_III/USUV_2025Jun24_NFLG_III_test.log') %>%
  pull(meanRate) %>%
  lapply(., AdjustFixedRate, 
         xml_filepath = "./2025Jun24/europe_clusters/workspace/USUV_2025Jun24_partial_III_SRD06_fixed_SG_updated.xml")




################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################