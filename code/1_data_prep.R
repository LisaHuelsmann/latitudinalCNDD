

# Prepare data and metadata




library(crayon)
library(tidyr)




# Define location of input and output -------------------------------------


## Input (raw data)
# if not stated otherwise, use example data in repo
# this is bci splitted up in half
if (!exists("path_input")) path_input = "data_prep/input/"


## Output (prepared data)
# if not stated otherwise, save in repo folders
if (!exists("path_output")) path_output = "data_prep/output/"




# Create output folders ---------------------------------------------------


# create folder for output
dir.create(path_output)

for (out in c("data_1_metainfo"
              , "data_3a_mortality"
              , "data_3b_growth"
              , "meta_abundances"
              , "meta_growth"
              , "meta_mortality"
              , "meta_stature")) {
  dir.create(paste0(path_output, out))
}




# Define sites ------------------------------------------------------------



# sites where tree or stem data is available
sites = unique(unlist(lapply(strsplit(c(list.files(paste0(path_input, "data_tree/"))
                                        , list.files(paste0(path_input, "data_stem/"))), 
                                      "_"), "[[", 1)))
               





# Prepare mortality data --------------------------------------------------


# Add metainfo
# Calculate competition
# Restructure for mortality

for (site in sites) {

  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  
  source("code/data_prep/data_1_metainfo.R")                        # produces 1_metainfo.Rdata

  source("code/data_prep/data_2_competition.R")                     # produces 2_competition.Rdata (not stored)

  source("code/data_prep/data_3a_mortality.R")                      # produces 3a_mortality.Rdata

}





# Prepare growth data -----------------------------------------------------


# Add metainfo
# no buffer applied
# Restructure for growth

for (site in sites) {
  
  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  
  source("code/data_prep/data_1_metainfo.R")                        # produces 1_metainfo.Rdata
  
  source("code/data_prep/data_3b_growth.R")                         # produces 3b_growth.Rdata
  
}




# Prepare meta ------------------------------------------------------------

# Calculates species-level summary statistics, that are saved in individual files 


# Abundances
for (site in sites) {
  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  source("code/data_prep/meta_abundances.R")
}


# Stature and maximum dbh
for (site in sites) {
  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  source("code/data_prep/meta_stature.R")
}


# Growth
for (site in sites) {
  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  source("code/data_prep/meta_growth.R")
}


# Mortality
for (site in sites) {
  flush.console()
  cat(blue(bold(paste(paste("\n", site, "\n")))))
  source("code/data_prep/meta_mortality.R")
}




