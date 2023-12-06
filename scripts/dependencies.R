###################################################################################################
###################################################################################################
# Dependencies for USUV European drivers analysis
###################################################################################################
############################# Check necessary packages are available and install; import if not.
pkgs <- c('tidyverse', 'magrittr' ,'ape', 'bioseq',  'parallel', 'ggpmisc', 'ggpubr', 'ggprism', 'cowplot',
         'scales','extrafont', 'future', 'rentrez', 'genbankr') 

for (pkg in pkgs){
  
  if (!(pkg %in% installed.packages())){
    
    install.packages(pkg)
    
    library(pkg, character.only = T)
    
  }else if (pkg %in% installed.packages()){
    
    library(pkg, character.only = T)
  }
}

cat('Dependencies Installed')
###################################################################################################
# Create directories for results and figures
yyyymonthdd <- format(Sys.Date(), '%Y%b%d')
check_dirs <- paste(c('./results', './figures'), yyyymonthdd, sep = '/')
dirs <- list.dirs()

for (check_dir in check_dirs){
  
  if (!(check_dir %in% dirs)){
    dir.create(check_dir, recursive = T)
  }
  
}

results_dir <- check_dirs[[1]]
figs_dir <- check_dirs[[2]]

###################################################################################################

# Execute in parallel
RunParallel <- function(func, v1, ...){
  options(warn = 1)
  require(parallel)
  
  # Set up cluster (fork)
  cl <- detectCores() %>% `-` (2) 
  
  start <- Sys.time()
  print(start)
  
  para <- mclapply(v1,
                   func, 
                   ...,
                   mc.cores = cl,
                   mc.set.seed = FALSE) #child process has the same initial random number generator (RNG) state as the current R session
  
  end <- Sys.time()
  elapsed <- end-start
  print(end)
  print(elapsed)
  
  return(para)
}


###################################################################################################

