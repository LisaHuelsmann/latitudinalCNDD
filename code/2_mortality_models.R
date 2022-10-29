


# Run mortality models

# Includes the following steps:
# 1) Grid search over different density definitions
# 2) Grid search summary to identify density definition for main analysis
# 3) Main mortality models (include more output)
# 4) Mortality models with randomized mortality/conspecific density


# NOTE: Grid search is time consuming!



library(tidyr)
library(dplyr)
library(parallel)
library(pbapply)




# Define location of input and output -------------------------------------

# if not stated otherwise, load and save from/in repo folders


# Data prep input
if (!exists("path_input")) path_input = "data_prep/input/"

# Data prep output
if (!exists("path_output")) path_output = "data_prep/output/"

# Analysis outputs
if (!exists("path_mortality")) path_mortality = "out/mortality_models/"


# create folder for output
dir.create(path_mortality)





# Define sites ------------------------------------------------------------


# sites where tree or stem data is available
sites = unique(unlist(lapply(strsplit(c(list.files(paste0(path_input, "data_tree/"))
                                        , list.files(paste0(path_input, "data_stem/"))), 
                                      "_"), "[[", 1)))





# Grid search mortality models --------------------------------------------


# define model runs
runs = c("gridsearchBAexpn", "gridsearchBAexp", "gridsearchNexpn", "gridsearchNexp")


for (run in runs) {
  
  # decay definitions
  decay_def = seq(1, 25, 2)
  decays_con = decays_tot = decay_def
  
  
  # chose n cores
  ncpu <- detectCores()-1
  ncpu <- ifelse(ncpu > 10, 3*length(decay_def), ncpu)
  
  
  # settings
  settings = expand.grid(decay_tot = decays_tot
                         , decay_con = decays_con
                         , site = sites, stringsAsFactors = F)
  settings = split(settings, seq(nrow(settings)))
  
  
  # create directory
  dir.create(paste0(path_mortality, run))
  for (decay_con in decays_con) {
    for (decay_tot in decays_tot) {
      dir.create(paste0(path_mortality, run, "/"
                        ,  sprintf("%02d", decay_con), "_", sprintf("%02d", decay_tot)))
    }
  }
  
  cl <- makeCluster(ncpu)
  source_setting = function(setting, run = run, paths) {
    site = setting$site
    decay_con = setting$decay_con
    decay_tot = setting$decay_tot
    source(paste0("code/mortality_models/mortality_", run, ".R"), local = T)
  }
  result <- pblapply(cl = cl, settings, source_setting, run = run,
                     paths = list(path_output = path_output, 
                                  path_input = path_input, 
                                  path_mortality = path_mortality))
  stopCluster(cl)
  
  
  # Combine results
  coefs_combined = data.frame()
  sums_combined = data.frame()
  
  
  for (decay_tot in decays_tot) {
    
    for (decay_con in decays_con) {
      
      # objects for global results
      coefs_global = data.frame()
      sums_global = data.frame()
      
      # vector of output object names
      output_objects = gsub("_global", "", ls()[grepl("_global", ls())])
      
      for (site in sites) {
        
        # load output objects
        load(paste0(path_mortality, run, "/"
                    , sprintf("%02d", decay_con), "_", sprintf("%02d", decay_tot), "/"
                    , site, "_mortality.Rdata"))
        
        # combine output objects
        for (i in output_objects) {
          temp = get(i)
          temp$site = site
          assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp)) 
        }
      }
      
      # combine different decay settings
      for (i in output_objects) {
        temp = get(paste0(i, "_global"))
        temp$decay_con = decay_con
        temp$decay_tot = decay_tot
        
        assign(paste0(i, "_combined"), rbind(get(paste0(i, "_combined")), temp))
      }
      
    }
    
  }
  
  
  # Save result 
  save(list = ls()[grepl("_combined", ls()) | ls() == "sites_ordered"]
       , file = paste0(path_mortality, run, "/combined_mortality.Rdata"))

}






# Grid search summary -----------------------------------------------------


# identifies optimal density definitions
# results stored in folder gridsearch_summary

source("code/mortality_models/grid_search_summary.R")






# Main mortality models ---------------------------------------------------


# optimal density definitions for these models are obtained from analyzing the grid search results
# see folder grid_search_summary


# chose model run
run = "main"

# chose n cores
ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 10, length(sites), ncpu)

# create directory
dir.create(paste0(path_mortality, run))


cl <- makeCluster(ncpu)
source_site = function(site, run, paths = paths) {
  source(paste0("code/mortality_models/mortality_", run, ".R"), local = T)
}

result <- pblapply(cl = cl, sites, source_site, run = run, 
                   paths = list(path_output = path_output, 
                                path_input = path_input, 
                                path_mortality = path_mortality))
stopCluster(cl)



# Session info
sink(paste0(path_mortality, run, "/", "sessionInfo.txt"))
sessionInfo()
sink()



# Combine results

# objects for global results
coefs_global = data.frame()
sums_global = data.frame()
AME_global = data.frame()
AMEsamples_global = data.frame()
rAME_global = data.frame()
rAMEsamples_global = data.frame()
nsp_global = data.frame()

# vector of output object names
output_objects = gsub("_global", "", ls()[grepl("_global", ls())])

for (site in sites) {
  
  # load output objects
  load(paste0(path_mortality, run, "/", site, "_mortality.Rdata"))
  
  # combine output objects
  for (i in output_objects) {
    temp = get(i)
    temp$site = site
    assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp)) 
  }
}


# Save result
save(list = ls()[grepl("_global", ls())]
     , file = paste0(path_mortality, run, "/global_mortality.Rdata"))


rm(list = ls()[!(grepl("site", ls()) | grepl("path", ls()))])






# Randomization mortality models ------------------------------------------


# define model runs
runs = c("randomizedMort", "randomizedConD")

# chose n cores
ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 10, length(sites), ncpu)


for (run in runs) {

  # create directory
  dir.create(paste0(path_mortality, run))
  
  cl <- makeCluster(ncpu)
  source_site = function(site, run, paths = paths) {
    source(paste0("code/mortality_models/mortality_", run, ".R"), local = T)
  }
  result <- pblapply(cl = cl, sites, source_site, run = run,
                     paths = list(path_output = path_output, 
                                  path_input = path_input, 
                                  path_mortality = path_mortality))
  stopCluster(cl)
  
  
  
  # Combine results
  
  # objects for global results
  coefs_global = data.frame()
  sums_global = data.frame()
  AME_global = data.frame()
  AMEsamples_global = data.frame()
  rAME_global = data.frame()
  rAMEsamples_global = data.frame()
  nsp_global = data.frame()
  
  # vector of output object names
  output_objects = gsub("_global", "", ls()[grepl("_global", ls())])
  
  for (site in sites) {
    
    # load output objects
    load(paste0(path_mortality, run, "/", site, "_mortality.Rdata"))
    
    # combine output objects
    for (i in output_objects) {
      temp = get(i)
      temp$site = site
      assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp)) 
    }
  }
  
  # Save result
  save(list = ls()[grepl("_global", ls())]
       , file = paste0(path_mortality, run, "/global_mortality.Rdata"))
  
}

rm(list = ls()[!(grepl("site", ls()) | grepl("path", ls()))])








