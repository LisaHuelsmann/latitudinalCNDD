


# Analyze global CNDD patterns

# Includes the following steps:
# Run analysis script for main and randomized mortality models


library(tidyr)
library(dplyr)
library(parallel)
library(pbapply)
library(metafor)
library(devtools)



# Define location of input and output -------------------------------------

# if not stated otherwise, load and save from/in repo folders


# Data prep input
if (!exists("path_input")) path_input = "data_prep/input/"

# Data prep output
if (!exists("path_output")) path_output = "data_prep/output/"

# Analysis outputs
if (!exists("path_mortality")) path_mortality = "out/mortality_models/"
if (!exists("path_meta")) path_meta = "out/meta_models/"


# create folder for output
dir.create(path_meta)




# Settings ----------------------------------------------------------------


# chose runs
runs = c("main", "randomizedMort", "randomizedConD")
run_names = c("Original dataset",
              "Randomized tree status",
              "Randomized conspecific densities")

# number cpu
ncpu = length(runs)
ncpu_meta = 15




# Run analysis ------------------------------------------------------------


cl <- makeCluster(ncpu)
analyze_run = function(run, ncpu_meta, paths) {
  
  # source analysis functions
  source("code/meta_models/functions_meta_models.R", local = T)
  
  # analyze results
  source("code/meta_models/source_meta_models.R", local = T) # must be local to use run object
  
  # Session info
  sink(paste0(paths$path_meta, "/", run,  "/sessionInfo.txt"))
  sessionInfo()
  sink()
  
}
result <- pblapply(cl = cl, runs, analyze_run, 
                   ncpu_meta = ncpu_meta,
                   paths = list(path_output = path_output, 
                                path_input = path_input, 
                                path_mortality = path_mortality,
                                path_meta = path_meta))
stopCluster(cl)







# Several run comparison figure -------------------------------------------


# latitudes for plotting
latitude_names = c("Tropical", "Subtropical", "Temperate")
latitudes = c(11.75, 29.25, 45)


# Load fitted models, save in list res
res = list()
for (run in runs) {
  load(paste0(path_meta, run, "/metamodels.Rdata"))
  res[[run]] = mod_list_red
}



# Fig 2 with three runs
dir.create(paste0(path_meta, "comparison"))
pdf(paste0(path_meta, "comparison/", "main_vs_randomized.pdf")
    , height = 7)

# loop through CNDD definitions
for (i in 1:length(res[[1]])) {
  
  # get results for each run of definition i
  res_i = lapply(res, "[[", i)
  
  # generate predictions for all runs
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1")
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0")

  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim))))
  
  # plotting
  layout(matrix(c(1, 1, 1, 2, 3, 4), ncol = 3, byrow = T))
  par(las = 1, oma = c(3, 4, 0, 0), mar = c(2, 2, 3, 1))
  
  # only latitude model
  plot_latitude(preds_lat, names = run_names
                , labelsx = 1, labelsy = 1
                , ylim = ylim
                , col_run = c("black", "dodgerblue3", "firebrick")
                , panel = "a")
  
  # latitude abundance model
  plot_abundance(preds_abund, names = run_names
                 , latitude_names = latitude_names
                 , labelsx = 2, labelsy = 1
                 , ylim = ylim
                 , col_run = c("black", "dodgerblue3", "firebrick")
                 , panel = "b")
  
}
dev.off()








# Get specific numbers for certain run ------------------------------------



# chose model runs
run = "main"

# Load fitted models
load(paste0(path_meta, run, "/metamodels.Rdata"))


# Choose CNDD definition
type = "rAME"
change = "equilibrium"
sel = lapply(mod_list_red, function(x) x$type) == type & lapply(mod_list_red, function(x) x$change) == change
res = mod_list_red[[which(sel)]]
backtrans = get(paste0("backtrans_", type))

# calculate predictions at two latitudes of interest (mean species CNDD model)
latitudes = c(11.75, 45)
newDat = cbind(transLat(latitudes, ref_lat = res$ref_lat))
predictions = predict(res$mod1, newmods = newDat, transf = backtrans)

# print values
cbind(latitude = latitudes,
      pred = round(100 * predictions$pred, 2),
      ci.lb = round(100 * predictions$ci.lb, 2),
      ci.ub = round(100 * predictions$ci.ub, 2))



# calculate predictions for rare and common species in the tropics (abundance-mediated CNDD model)
latitudes = 11.75
abundances = c(1, 100)
newDat = cbind(transLat(latitudes, ref_lat = res$ref_lat)
               , transAbund(abundances, ref_abund = res$ref_abund)
               , transLat(latitudes, ref_lat = res$ref_lat)*transAbund(abundances, ref_abund = res$ref_abund))
predictions = predict(res$mod0, newmods = newDat, transf = backtrans)

# print values
cbind(abundances = abundances,
      latitudes = latitudes,
      pred = round(100 * predictions$pred, 2),
      ci.lb = round(100 * predictions$ci.lb, 2),
      ci.ub = round(100 * predictions$ci.ub, 2))


# calculate predictions for rare and common species in temperate forests (abundance-mediated CNDD model)
latitudes = 45
abundances = c(1, 100)
newDat = cbind(transLat(latitudes, ref_lat = res$ref_lat)
               , transAbund(abundances, ref_abund = res$ref_abund)
               , transLat(latitudes, ref_lat = res$ref_lat)*transAbund(abundances, ref_abund = res$ref_abund))
predictions = predict(res$mod0, newmods = newDat, transf = backtrans)

# print values
cbind(abundances = abundances,
      latitudes = latitudes,
      pred = round(100 * predictions$pred, 2),
      ci.lb = round(100 * predictions$ci.lb, 2),
      ci.ub = round(100 * predictions$ci.ub, 2))




# # calculate predictions for rare species in the tropics and temperate forests (abundance-mediated CNDD model)
# latitudes = c(11.75, 45)
# abundances = 1
# newDat = cbind(transLat(latitudes, ref_lat = res$ref_lat)
#                , transAbund(abundances, ref_abund = res$ref_abund)
#                , transLat(latitudes, ref_lat = res$ref_lat)*transAbund(abundances, ref_abund = res$ref_abund))
# predictions = predict(res$mod0, newmods = newDat, transf = backtrans)
# 
# # print values
# cbind(abundances = abundances,
#       latitudes = latitudes,
#       pred = round(100 * predictions$pred, 2),
#       ci.lb = round(100 * predictions$ci.lb, 2),
#       ci.ub = round(100 * predictions$ci.ub, 2))
# 
# 
# # calculate predictions for common species in the tropics and temperate forests (abundance-mediated CNDD model)
# latitudes = c(11.75, 45)
# abundances = 100
# newDat = cbind(transLat(latitudes, ref_lat = res$ref_lat)
#                , transAbund(abundances, ref_abund = res$ref_abund)
#                , transLat(latitudes, ref_lat = res$ref_lat)*transAbund(abundances, ref_abund = res$ref_abund))
# predictions = predict(res$mod0, newmods = newDat, transf = backtrans)
# 
# # print values
# cbind(abundances = abundances,
#       latitudes = latitudes,
#       pred = round(100 * predictions$pred, 2),
#       ci.lb = round(100 * predictions$ci.lb, 2),
#       ci.ub = round(100 * predictions$ci.ub, 2))


