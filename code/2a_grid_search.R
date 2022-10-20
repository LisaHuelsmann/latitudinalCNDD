

# comparison of decay values
# full factorial 


library(dplyr)   # for data wrangling
library(viridis) # for colors
library(fields)  # for image.plot



# Load data ---------------------------------------------------------------


# chose model runs
# NOTE:these must be based on the same data and particularly the same max distance (i.e. check columns Data_from and Max_thresh)


### for full exploration
settings = data.frame(run = c( "el", "ek", "ej", "ei")
                      , decay_type = c("exp", "expn", "exp", "expn")
                      , term = c("N", "N", "BA", "BA"))

# ### for subset
# settings = data.frame(run = c("ei", "ej")
#                       , decay_type = c("expn", "exp")
#                       , term = c("BA", "BA"))

# # for full and subset
# settings = data.frame(run = c("ec", "ed", "ee", "ef", "eg", "eh")
#                       , decay_type = c("expn", "expn", "exp", "exp", "expn", "exp")
#                       , term = c("N", "BA", "N", "BA", "BA", "BA"))


# object for all results
sums_all = data.frame()


for (s in 1:nrow(settings)) {
  
  # load combined data
  load(paste0("out/mortality_models/", settings$run[s], "/combined_mortality.Rdata"))
  temp = sums_combined
  
  # chose columns (not all outputs have the same)
  temp = temp[, c("logLik", "AIC", "nobs", "sp", "AUC", "site", "latitude", "abundance", "decay_con", "decay_tot")]
  
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
# YES, 169 combinations, 2500 species




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


pdf(paste0("out/mortality_models/", settings$run[1], "/sumLL_", paste(settings$run, collapse = "_"), ".pdf"), height = 9, width = 9)
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

# copy figure to all other run folders
for (i in 2:nrow(settings)) {
  file.copy(paste0("out/mortality_models/", settings$run[1], "/sumLL_", paste(settings$run, collapse = "_"), ".pdf")
            , paste0("out/mortality_models/", settings$run[i], "/")
            , overwrite = T)
}





# Plot mean of AUC --------------------------------------------------------


pdf(paste0("out/mortality_models/", settings$run[1], "/meanAUC_", paste(settings$run, collapse = "_"), ".pdf"), height = 9, width = 9)
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

# copy figure to all other run folders
for (i in 2:nrow(settings)) {
  file.copy(paste0("out/mortality_models/", settings$run[1], "/meanAUC_", paste(settings$run, collapse = "_"), ".pdf")
            , paste0("out/mortality_models/", settings$run[i], "/")
            , overwrite = T)
}



# Save global optimum -----------------------------------------------------


sink(paste0("out/mortality_models/", settings$run[1], "/globalOpt_", paste(settings$run, collapse = "_"), ".txt"))
for (criterion in c("sumlogLik", "meanAUC")) {
  opt = sums_total[which.max(sums_total[, criterion]), ]
  cat(criterion, "is optimal at", "\n")
  cat("term =", opt$term, "\n")
  cat(opt$decay_type, "decay", "\n")
  cat(paste(c("conspecific", "total"), "decay =", opt[, c("decay_con", "decay_tot")], collapse = "\n"), "\n", "\n")
}
sink()

# copy figure to all other run folders
for (i in 2:nrow(settings)) {
  file.copy(paste0("out/mortality_models/", settings$run[1], "/globalOpt_", paste(settings$run, collapse = "_"), ".txt")
            , paste0("out/mortality_models/", settings$run[i], "/")
            , overwrite = T)
}





# Compare AUC at AUC and LL optimum ---------------------------------------


opt_AUC = sums_total[which.max(sums_total[, "meanAUC"]), ]$meanAUC
opt_LL = sums_total[which.max(sums_total[, "sumlogLik"]), ]$meanAUC
opt_AUC-opt_LL



# Site specific optimum ---------------------------------------------------


sums_all_mod %>% 
  filter(!site_sp %in% incomplete) %>%
  group_by(site, latitude, run, term, decay_con, decay_tot, decay_type) %>% 
  summarise(sumlogLik = sum(logLik),
            meanlogLik = mean(logLik), 
            meanAUC = mean(AUC, na.rm = T),
            medianAUC = median(AUC, na.rm = T),
            nvalues = n()) %>% 
  group_by(site) -> sums_total_site


sums_total_site %>% 
  filter(sumlogLik == max(sumlogLik)) -> site_opt

# there is quite some variability between sites
hist(site_opt$decay_con[site_opt$decay_type == "exp"])
hist(site_opt$decay_con[site_opt$decay_type == "expn"])
hist(site_opt$decay_tot[site_opt$decay_type == "exp"])
hist(site_opt$decay_tot[site_opt$decay_type == "expn"])
plot(decay_con ~ nvalues, site_opt)
plot(decay_con ~ decay_tot, site_opt)
plot(decay_con ~ latitude, site_opt)


# there is no pattern with latitude



# Species specific optimum ---------------------------------------------------

sums_all_mod %>% 
  filter(!site_sp %in% incomplete) %>%
  group_by(site_sp, decay_type) %>% 
  filter(logLik == max(logLik)) -> sums_total_species

opt_species_exp = sums_total_species[sums_total_species$decay_type == "exp", ]
opt_species_expn = sums_total_species[sums_total_species$decay_type == "expn", ]

smoothScatter(opt_species_exp$decay_con, opt_species_exp$decay_tot)
hist(opt_species_exp$decay_con)
hist(opt_species_expn$decay_con)

scatter.smooth(opt_species_exp$decay_con, opt_species_exp$logLik)
scatter.smooth(opt_species_expn$decay_con, opt_species_expn$logLik)
scatter.smooth(opt_species_exp$decay_con, opt_species_exp$logLik, ylim = c(-1000, 0))


scatter.smooth(log(opt_species_exp$abundance), opt_species_exp$logLik)
scatter.smooth(log(opt_species_expn$abundance), opt_species_expn$logLik)

