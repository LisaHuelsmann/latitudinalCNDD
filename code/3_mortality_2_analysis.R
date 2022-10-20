


# Analyze CNDD in mortality
# sources analysis scripts


library(tidyr)
library(dplyr)
library(parallel)
library(pbapply)
library(metafor)
library(devtools)



# Settings ----------------------------------------------------------------


# Define outstyle
outstyle = "j"

# chose model runs
# runs = c("ax", "fa", "fb", "fc")
# runs = c("zb", "zc")
# runs = c("zc", "fd", "fe", "ff")
# runs = c("zd")
# runs = c("fg", "fh")
# runs = "zf"
runs = c("fi")
# runs = c("fj")
# runs = c("zf", "fi", "fj")

# numper cpu
ncpu = length(runs)
ncpu_meta = 6




# Run analysis ------------------------------------------------------------


cl <- makeCluster(ncpu)
analyze_run = function(run, outstyle, ncpu_meta) {
  
  # turbo colors
  devtools::source_gist("https://gist.github.com/jlmelville/be981e2f36485d8ef9616aef60fd52ab")
  
  # source analysis functions
  source("code/functions_analysis.R", local = T)
  
  # analyze results
  source(paste0("code/__mortality_2", outstyle, "_analysis.R"), local = T) # must be local to use run object
  
  # Session info
  sink(paste0("out/mortality_models/", run,  "/analysis_", outstyle, "/", "sessionInfo.txt"))
  sessionInfo()
  sink()
  
}
result <- pblapply(cl = cl, runs, analyze_run, outstyle = outstyle, ncpu_meta = ncpu_meta)
stopCluster(cl)



# for h analysis and ay run:         51m 17s (I think on DELL?)
# for i analysis and ay run:     01h 04m 16s (on DELL)
# for h analysis and az run:         43m 54s (on server)
# for h analysis and za run:         44m 03s (on server)
# for h analysis and zb run:         17m 06s (on MAC) ncpu_meta = 3
# for h analysis and 4 runs:         24m 09s (on UBT server) ncpu_meta = 3
# for h analysis and zc runs:        24m 09s (on UBT server) ncpu_meta = 3
# for h analysis and fe ff runs:     19m 47s (on UBT server) ncpu_meta = 8
# for j analysis and zf runs:    04h 14m 17s (on UBT server) ncpu_meta = 6











# Several run comparison figure -------------------------------------------



# Define outstyle
outstyle = "j" # must be the same for all because of transformation functions

# chose model runs
runs = c("zf", "fi", "fj")

# names for runs
run_names = c("Original dataset",
              "Randomized conspecific densities",
              "Randomized tree status"
              )

# latitudes for plotting
latitude_names = c("Tropical", "Subtropical", "Temperate")
latitudes = c(11.75, 29.25, 45)


# Load fitted models, save in list res
res = list()
for (run in runs) {
  load(paste0("out/mortality_models/", run, "/analysis_", outstyle, "/metamodels.Rdata"))
  res[[run]] = mod_list_red
}



# Fig 2 with three runs
pdf(paste0("out/mortality_models/comparison_", paste(runs, collapse = "_"), "_outstyle_", outstyle, ".pdf")
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
  par(las = 1, oma = c(3, 4, 3, 0), mar = c(2, 2, 3, 1))
  
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



# Define outstyle
outstyle = "j"

# chose model runs
run = "zf"

# Load fitted models
load(paste0("out/mortality_models/", run, "/analysis_", outstyle, "/metamodels.Rdata"))


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
abundances = c(1, 1000)
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
abundances = c(1, 1000)
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




# calculate predictions for rare species in the tropics and temperate forests (abundance-mediated CNDD model)
latitudes = c(11.75, 45)
abundances = 1
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


# calculate predictions for common species in the tropics and temperate forests (abundance-mediated CNDD model)
latitudes = c(11.75, 45)
abundances = 1000
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


