

# Grid search comparison to identify optimal density definitions


library(dplyr)   # for data wrangling
library(viridis) # for colors


# create directory
dir.create(paste0(path_mortality, "gridsearch_summary"))



# Load data ---------------------------------------------------------------


# chose model runs
# NOTE:these must be based on the same data and particularly the same max distance


### for full exploration
settings = data.frame(runs = c("gridsearchNexp",  "gridsearchNexpn", "gridsearchBAexp", "gridsearchBAexpn")
                      , decay_type = c("exp", "expn", "exp", "expn")
                      , term = c("N", "N", "BA", "BA"))


# object for all results
sums_all = data.frame()


for (s in 1:nrow(settings)) {
  
  # load combined data
  load(paste0(path_mortality, settings$run[s], "/combined_mortality.Rdata"))
  temp = sums_combined
  
  # chose columns (not all outputs have the same)
  temp = temp[, c("logLik", "AIC", "nobs", "sp", "AUC", "site", "decay_con", "decay_tot")]
  
  # add term, run and decay type
  temp$run = settings$run[s]
  temp$term = settings$term[s]
  temp$decay_type = settings$decay_type[s]
  
  sums_all = rbind(sums_all, temp)
  
  rm(list = c("temp", ls()[grepl("_combined", ls())]))
  
}

# explore available results
table(sums_all$run)
table(sums_all$site, sums_all$term)
table(sums_all$site, sums_all$decay_type)

table(sums_all$decay_tot, sums_all$term, sums_all$decay_type)
table(sums_all$decay_con, sums_all$term, sums_all$decay_type)





# Prepare data ------------------------------------------------------------


# add site_sp column
sums_all %>% 
  mutate(site_sp = paste(site, sp, sep = "_"),
         run_site_sp = paste(run, site, sp, sep = "_"),
         decay_combo = paste(decay_con, decay_tot, sep = "_")) -> sums_all_mod

# get incomplete run-site-species combinations
sums_all_mod %>% 
  group_by(run_site_sp) %>% 
  mutate(ncombos = length(unique(decay_combo))) %>%              # how many combinations of decay values per run
  group_by(run) %>%           
  mutate(maxcombos = max(ncombos)) %>%                           # how many combinations of decay values maximum per run
  group_by(site_sp) %>%       
  summarize(run_problem = length(unique(run)) < nrow(settings),  # less runs then expected per site_sp
            combo_problem = ncombos < maxcombos) -> counts       # less decay combos than expected

incomplete = counts$site_sp[counts$run_problem | counts$combo_problem]




# Aggregate different criteria per setting --------------------------------


sums_all_mod %>% 
  filter(!site_sp %in% incomplete) %>%
  group_by(run, term, decay_con, decay_tot, decay_type) %>%
  summarise(nvalues = n(),
            sumlogLik = sum(logLik),
            meanlogLik = mean(logLik), 
            meanAUC = mean(AUC, na.rm = T),
            medianAUC = median(AUC, na.rm = T)) %>% 
  arrange(run, decay_con, decay_tot) %>%
  as.data.frame() -> sums_total

# check if really the same number of models was used
table(sums_total$run, sums_total$nvalues)
# YES, 169 combinations




# Functions ---------------------------------------------------------------


# plot against setting
plot_map = function(run, term, type, criterion, optimum = T) {
  
  # prep matrix
  dat = sums_total[sums_total$run == run & sums_total$term == term & sums_total$decay_type == type, ]
  mat = matrix(dat[, criterion], ncol = length(unique(dat$decay_tot)), byrow = T)
  
  
  # plot matrix
  image(unique(dat$decay_con), unique(dat$decay_tot), mat
        , col = viridis(128)
        , zlim = range(sums_total[, criterion])
        , main = paste0(term, ", ", ifelse(type == "exp", "exponential", "exponential-normal"), " decay")
        , axes = F
        , asp = 1
        , xlab = ""
        , ylab = ""
        , bty='l'
  )
  axis(1, at = seq(0, 25, 5), pos = 0)
  axis(2, at = seq(0, 25, 5), pos = 0)
  mtext(expression(paste(mu, " con")), side = 1, line = 2.1)
  mtext(expression(paste(mu, " tot")), side = 2, line = 2, las = 0)
  
  # add optimum
  if (optimum) {
    opt = dat[which.max(dat[, criterion]), ]
    points(opt$decay_con, opt$decay_tot, col = "black", pch = 16)
    text(opt$decay_con, opt$decay_tot-1.5, "opt", col = "black", cex = 0.8)
  }
}




# Plot sum of log likelihood ----------------------------------------------


pdf(paste0(path_mortality, "/gridsearch_summary/sumLL.pdf"), height = 9, width = 9)
par(mfrow = c(ceiling(sqrt(nrow(settings))), floor(sqrt(nrow(settings)))), mar = c(4, 4, 2, 1), las = 1)

# get global optimum
opt = sums_total[which.max(sums_total[, "sumlogLik"]), ]

for (s in 1:nrow(settings)) {
  plot_map(run = settings$run[s]
           , term = settings$term[s]
           , type = settings$decay_type[s]
           , criterion = "sumlogLik"
           , optimum = settings$run[s] == opt$run)
}
dev.off()




# Plot mean of AUC --------------------------------------------------------


pdf(paste0(path_mortality, "/gridsearch_summary/meanAUC.pdf"), height = 9, width = 9)
par(mfrow = c(ceiling(sqrt(nrow(settings))), floor(sqrt(nrow(settings)))), mar = c(4, 4, 2, 1), las = 1)

# get global optimum
opt = sums_total[which.max(sums_total[, "meanAUC"]), ]

for (s in 1:nrow(settings)) {
  plot_map(run = settings$run[s]
           , term = settings$term[s]
           , type = settings$decay_type[s]
           , criterion = "meanAUC"
           , optimum = settings$run[s] == opt$run)
}
dev.off()



# Save global optimum -----------------------------------------------------


sink(paste0(path_mortality, "/gridsearch_summary/globalOpt.txt"))
for (criterion in c("sumlogLik", "meanAUC")) {
  opt = sums_total[which.max(sums_total[, criterion]), ]
  cat(criterion, "is optimal at", "\n")
  cat("term =", opt$term, "\n")
  cat(opt$decay_type, "decay", "\n")
  cat(paste(c("conspecific", "total"), "decay =", opt[, c("decay_con", "decay_tot")], collapse = "\n"), "\n", "\n")
}

# AUC difference at AUC and LL optimum
opt_AUC = sums_total[which.max(sums_total[, "meanAUC"]), ]$meanAUC
opt_LL = sums_total[which.max(sums_total[, "sumlogLik"]), ]$meanAUC
cat("AUC difference between AUC and LL optimum =", round(opt_AUC-opt_LL, 3), "\n") 

# Number of species considered
cat("number of species =", unique(sums_total$nvalues))

sink()


# Clean

rm(list = ls()[!(grepl("site", ls()) | grepl("path", ls()))])





