


# Run models for CNDD in mortality



library(tidyr)
library(dplyr)
library(parallel)
library(pbapply)






# Define sites ------------------------------------------------------------

# loop through all sites where tree or stem data is available
sites = sort(unique(unlist(lapply(strsplit(c(list.files("../ForestGEO_datacleaning@git/data_tree/")
                                             , list.files("../ForestGEO_datacleaning@git/data_stem/")), 
                                           "_"), "[[", 1)))
             , decreasing = F)

# subset for trials
# sites = sites[10:15]
# sites = "serc"





# Type a and z ------------------------------------------------------------

# Run models --------------------------------------------------------------



# chose model run
run = "zf"


dir.create(paste0("out/mortality_models/", run))

ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 4, 20, ncpu)
cl <- makeCluster(ncpu)
source_site = function(site, run = run) {
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, sites, source_site, run = run)
stopCluster(cl)



# Session info
sink(paste0("out/mortality_models/", run, "/", "sessionInfo.txt"))
sessionInfo()
sink()





# Load and prepare results ------------------------------------------------


# objects for global results
coefs_global = data.frame()
sums_global = data.frame()
AME_global = data.frame()
AMEsamples_global = data.frame()
rAME_global = data.frame()
rAMEsamples_global = data.frame()
nsp_global = data.frame()


for (site in sites) {
  
  # load output objects
  load(paste0("out/mortality_models/", run, "/", site, "_mortality.Rdata"))
  
  # vector of output object names
  output_objects = gsub("_global", "", ls()[grepl("_global", ls())])
  
  # add meta information (species and site level)
  source("code/source_add_meta_information.R")
  
  rm(list = c("temp", output_objects))
}


# object with sites ordered by absolute latitude
AME_global %>% 
  group_by(site) %>% 
  summarise(latitude = unique(latitude)) %>% 
  arrange(abs(latitude)) -> sites_ordered





# Save result -------------------------------------------------------------


save(list = ls()[grepl("_global", ls()) | ls() == "sites_ordered"]
     , file = paste0("out/mortality_models/", run, "/global_mortality.Rdata"))










# Type b ------------------------------------------------------------------


# Run models --------------------------------------------------------------


# chose model run
run = "bb"

# exponent
exponent = seq(0.01, 1.4, length.out = 11)

# settings
settings = expand.grid(exponent = exponent, site = sites, stringsAsFactors = F)
settings = split(settings, seq(nrow(settings)))

# create directory
dir.create(paste0("out/mortality_models/", run))
for (e in exponent) dir.create(paste0("out/mortality_models/", run, "/", format(e, nsmall = 3)))

ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 15, length(exponent), ncpu)
cl <- makeCluster(ncpu)
source_setting = function(setting, run = run) {
  site = setting$site
  e = setting$exponent
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, settings, source_setting, run = run)
stopCluster(cl)



# Load and prepare results ------------------------------------------------

# table with meta data
Site_table = as.data.frame(readxl::read_excel("../ForestGEO_datacleaning@git/global_metainfo/Sites.xlsx", sheet = 1, skip = 1))


for (e in exponent) {
  # objects for global results
  coefs_global = data.frame()
  sums_global = data.frame()
  AME_global = data.frame()
  RR_global = data.frame()
  nsp_global = data.frame()
  
  for (site in sites) {
    
    # load output objects
    load(paste0("out/mortality_models/", run, "/", format(e, nsmall = 3), "/", site, "_mortality.Rdata"))
    load(paste0("data_prep/meta_abundances/", site, "_abundances_tree.Rdata"))
    
    # add additional info to all output style objects
    for (i in gsub("_global", "", ls()[grepl("_global", ls())])) {
      
      temp = get(i)
      temp$site = site
      
      # add latitude
      temp$latitude = Site_table$Latitude[match(toupper(temp$site), Site_table$Abbreviation)]
      
      # add abundances per species
      temp$abundance = abundances_tree$Nha[match(temp$sp, abundances_tree$sp)]  
      
      # add mean abundance for rare species
      temp$abundance[temp$sp == "Rare_tree"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                          nsp$sp[which(nsp$rare_stature == "Rare_tree")]])
      temp$abundance[temp$sp == "Rare_shrub"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                           nsp$sp[which(nsp$rare_stature == "Rare_shrub")]])
      
      assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp))
    }
    rm(list = c("temp", gsub("_global", "", ls()[grepl("_global", ls())])))
  }
  
  
  # object with sites ordered by absolute latitude
  AME_global %>% 
    group_by(site) %>% 
    summarise(latitude = unique(latitude)) %>% 
    arrange(abs(latitude)) -> sites_ordered
  
  
  # remove SE = 0 or NA
  sel = AME_global$std.error !=0 & !is.na(AME_global$std.error)
  AME_global = AME_global[sel, ]
  RR_global = RR_global[sel, ]
  
  
  
  
  # Save result -------------------------------------------------------------
  
  
  save(list = ls()[grepl("_global", ls()) | ls() == "sites_ordered"]
       , file = paste0("out/mortality_models/", run, "/", format(e, nsmall = 3), "/global_mortality.Rdata"))

}


rm(list = ls())




# Type c ------------------------------------------------------------------

# Run models --------------------------------------------------------------


# chose model run
run = "cg"

# distance definitions
dist_def = seq(4, 30, 2)

# settings
settings = expand.grid(dist_def = dist_def, site = sites, stringsAsFactors = F)
settings = split(settings, seq(nrow(settings)))

# create directory
dir.create(paste0("out/mortality_models/", run))
for (distance in dist_def) dir.create(paste0("out/mortality_models/", run, "/", sprintf("%02d", distance)))


ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 15, length(dist_def), ncpu)
cl <- makeCluster(ncpu)
source_setting = function(setting, run = run) {
  site = setting$site
  distance = setting$dist_def
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, settings, source_setting, run = run)
stopCluster(cl)







# Load and prepare results ------------------------------------------------

# table with meta data
Site_table = as.data.frame(readxl::read_excel("../ForestGEO_datacleaning@git/global_metainfo/Sites.xlsx", sheet = 1, skip = 1))


for (distance in dist_def) {
  
  # objects for global results
  coefs_global = data.frame()
  sums_global = data.frame()
  nsp_global = data.frame()
  AME_global = data.frame()
  rAME_global = data.frame()
  
  for (site in sites) {
    
    # load output objects
    load(paste0("out/mortality_models/", run, "/", sprintf("%02d", distance), "/", site, "_mortality.Rdata"))
    load(paste0("data_prep/meta_abundances/", site, "_abundances_tree.Rdata"))
    
    # add additional info to all output style objects
    for (i in gsub("_global", "", ls()[grepl("_global", ls())])) {
      
      temp = get(i)
      temp$site = site
      
      # add latitude
      temp$latitude = Site_table$Latitude[match(toupper(temp$site), Site_table$Abbreviation)]
      
      # add abundances per species
      temp$abundance = abundances_tree$Nha[match(temp$sp, abundances_tree$sp)]  
      
      # add mean abundance for rare species
      temp$abundance[temp$sp == "Rare_tree"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                          nsp$sp[which(nsp$rare_stature == "Rare_tree")]])
      temp$abundance[temp$sp == "Rare_shrub"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                           nsp$sp[which(nsp$rare_stature == "Rare_shrub")]])
      
      assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp))
    }
    rm(list = c("temp", gsub("_global", "", ls()[grepl("_global", ls())])))
  }
  
  
  
  # Save result -------------------------------------------------------------
  
  
  save(list = ls()[grepl("_global", ls()) | ls() == "sites_ordered"]
       , file = paste0("out/mortality_models/", run, "/", sprintf("%02d", distance), "/global_mortality.Rdata"))
  
}







# Type d ------------------------------------------------------------------

# Run models --------------------------------------------------------------


# chose model run
run = "db"

# decay definitions
decay_def = seq(0, 1, 0.1)

# settings
settings = expand.grid(decay_def = decay_def, site = sites, stringsAsFactors = F)
settings = split(settings, seq(nrow(settings)))

# create directory
dir.create(paste0("out/mortality_models/", run))
for (decay in decay_def) dir.create(paste0("out/mortality_models/", run, "/", format(decay, nsmall = 1)))


ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 15, length(decay_def), ncpu)
cl <- makeCluster(ncpu)
source_setting = function(setting, run = run) {
  site = setting$site
  decay = setting$decay_def
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, settings, source_setting, run = run)
stopCluster(cl)







# Load and prepare results ------------------------------------------------

# table with meta data
Site_table = as.data.frame(readxl::read_excel("../ForestGEO_datacleaning@git/global_metainfo/Sites.xlsx", sheet = 1, skip = 1))


for (decay in decay_def) {
  
  # objects for global results
  coefs_global = data.frame()
  sums_global = data.frame()
  nsp_global = data.frame()
  AME_global = data.frame()
  rAME_global = data.frame()
  
  for (site in sites) {
    
    # load output objects
    load(paste0("out/mortality_models/", run, "/", format(decay, nsmall = 1), "/", site, "_mortality.Rdata"))
    load(paste0("data_prep/meta_abundances/", site, "_abundances_tree.Rdata"))
    
    # add additional info to all output style objects
    for (i in gsub("_global", "", ls()[grepl("_global", ls())])) {
      
      temp = get(i)
      temp$site = site
      
      # add latitude
      temp$latitude = Site_table$Latitude[match(toupper(temp$site), Site_table$Abbreviation)]
      
      # add abundances per species
      temp$abundance = abundances_tree$Nha[match(temp$sp, abundances_tree$sp)]  
      
      # add mean abundance for rare species
      temp$abundance[temp$sp == "Rare_tree"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                          nsp$sp[which(nsp$rare_stature == "Rare_tree")]])
      temp$abundance[temp$sp == "Rare_shrub"] = mean(abundances_tree$Nha[abundances_tree$sp %in% 
                                                                           nsp$sp[which(nsp$rare_stature == "Rare_shrub")]])
      
      assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp))
    }
    rm(list = c("temp", gsub("_global", "", ls()[grepl("_global", ls())])))
  }
  
  
  
  # Save result -------------------------------------------------------------
  
  
  save(list = ls()[grepl("_global", ls()) | ls() == "sites_ordered"]
       , file = paste0("out/mortality_models/", run, "/", format(decay, nsmall = 1), "/global_mortality.Rdata"))
  
}





# Type e ------------------------------------------------------------------

# Run models --------------------------------------------------------------


# chose model run
run = "el"

# decay definitions
decay_def = seq(1, 25, 2)
decays_con = decays_tot = decay_def
# decays_con = seq(1, 4.5, 0.5)
# decays_tot = seq(17, 24, 1)



# settings
settings = expand.grid(decay_tot = decays_tot
                       , decay_con = decays_con
                       , site = sites, stringsAsFactors = F)
settings = split(settings, seq(nrow(settings)))


# create directory
dir.create(paste0("out/mortality_models/", run))
for (decay_con in decays_con) {
  for (decay_tot in decays_tot) {
  dir.create(paste0("out/mortality_models/", run, "/"
                    ,  format(decay_con, nsmall = 1), "_", sprintf("%02d", decay_tot)))
  }
}


ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 15, 26, ncpu)
cl <- makeCluster(ncpu)
source_setting = function(setting, run = run) {
  site = setting$site
  decay_con = setting$decay_con
  decay_tot = setting$decay_tot
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, settings, source_setting, run = run)
stopCluster(cl)






# Load and prepare results ------------------------------------------------


coefs_combined = data.frame()
sums_combined = data.frame()


for (decay_tot in decays_tot) {
  
  for (decay_con in decays_con) {
    
    # objects for global results
    coefs_global = data.frame()
    sums_global = data.frame()
    
    for (site in sites) {
      
      # load output objects
      load(paste0("out/mortality_models/", run, "/"
                  , format(decay_con, nsmall = 1), "_", sprintf("%02d", decay_tot), "/"
                  , site, "_mortality.Rdata"))

      # vector of output object names
      output_objects = gsub("_global", "", ls()[grepl("_global", ls())])
      
      # add meta information (species and site level)
      source("code/source_add_meta_information.R")
      
      rm(list = c("temp", output_objects))
      
    }
    
    for (i in output_objects) {
      temp = get(paste0(i, "_global"))
      temp$decay_con = decay_con
      temp$decay_tot = decay_tot
      
      assign(paste0(i, "_combined"), rbind(get(paste0(i, "_combined")), temp))
    }

  }
  
}


# Save result -------------------------------------------------------------


save(list = ls()[grepl("_combined", ls()) | ls() == "sites_ordered"]
     , file = paste0("out/mortality_models/"
                     , run, "/combined_mortality.Rdata"))





# Type f ------------------------------------------------------------------


# Run models --------------------------------------------------------------



# chose model run
run = "fj"


dir.create(paste0("out/mortality_models/", run))

ncpu <- detectCores()-1
ncpu <- ifelse(ncpu > 20, 12, ncpu)
cl <- makeCluster(ncpu)
source_site = function(site, run = run) {
  source(paste0("code/mortality_", run, ".R"), local = T)
}
result <- pblapply(cl = cl, sites, source_site, run = run)
stopCluster(cl)




# Load and prepare results ------------------------------------------------


# objects for global results
coefs_global = data.frame()
sums_global = data.frame()
AME_global = data.frame()
AMEsamples_global = data.frame()
rAME_global = data.frame()
rAMEsamples_global = data.frame()
nsp_global = data.frame()

for (site in sites) {
  
  # load output objects
  load(paste0("out/mortality_models/", run, "/", site, "_mortality.Rdata"))

  # vector of output object names
  output_objects = gsub("_global", "", ls()[grepl("_global", ls())])
  
  # add meta information (species and site level)
  source("code/source_add_meta_information.R")

  rm(list = c("temp", output_objects))
}


# object with sites ordered by absolute latitude
AME_global %>% 
  group_by(site) %>% 
  summarise(latitude = unique(latitude)) %>% 
  arrange(abs(latitude)) -> sites_ordered





# Save result -------------------------------------------------------------


save(list = ls()[grepl("_global", ls()) | ls() == "sites_ordered"]
     , file = paste0("out/mortality_models/", run, "/global_mortality.Rdata"))


















# Cleaning outputs --------------------------------------------------------


# remove unnessary files
# Rdata files for individual sites (for all runs)
# pdfs with gam splines (from b runs on)






# # xxx_mortality.Rdata can be deleted if global_mortality.Rdata exists
# runs = list.files("out/mortality_models/")
# runs = runs[!grepl("\\.", runs)]
# 
# for (i in runs) {
# 
#   path = paste0("out/mortality_models/", i)
#   content = list.files(path, full.names = T) # all content
# 
#   # some of the content has global_mortality.R in name --> results are in main folder
#   if (any(grepl("global_mortality.R", content))) {  # deletion procedure in main folder
#     paths = path
#   } else {                                          # deletion procedure in subfolders
#     paths = list.dirs(path, recursive = F)
#   }
# 
#   for (j in paths) {
# 
#     # list Rdata files
#     files = list.files(j)
#     files = files[grepl(".Rdata", files)]
# 
#     # does global_mortality.Rdata exist?
#     if (any(files == "global_mortality.Rdata")) {
# 
#       # delete files except for global_mortality.Rdata
#       sapply(paste0(j, "/", files[!grepl("global_mortality.Rdata", files)]), unlink)
#     }
# 
#     # if not, no deletion
# 
#   }
# 
# }








# # xxx.pdf for grid search runs (b-e) take too much space
# runs = list.files("out/mortality_models/")
# runs = runs[!grepl("\\.", runs)]
# runs = runs[!startsWith(runs, "a") & !startsWith(runs, "z")]
# 
# for (i in runs) {
# 
#   path = paste0("out/mortality_models/", i)
#   content = list.files(path, full.names = T) # all content
# 
#   # some of the content has global_mortality.R in name --> results are in main folder
#   if (any(grepl("global_mortality.R", content))) {  # deletion procedure in main folder
#     paths = path
#   } else {                                          # deletion procedure in subfolders
#     paths = list.dirs(path, recursive = F)
#   }
# 
#   for (j in paths) {
# 
#     # list pdf files
#     files = list.files(j)
#     files = files[grepl(".pdf", files)]
# 
#     # delete files
#     sapply(paste0(j, "/", files), unlink)
# 
#   }
# 
# }





