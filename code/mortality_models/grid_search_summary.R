

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


# Meta info sites
Site_table = readxl::read_xlsx(paste0(path_input, "global_metainfo/plot_sites_information.xlsx"), sheet = 1)




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




# Species-specific optima -------------------------------------------------


sums_all_mod %>% 
  filter(!site_sp %in% incomplete) %>%
  filter(run == opt$run) %>%  # only run that was optimal overall (i.e BA or N, exp or expn)
  group_by(site_sp) %>% 
  arrange(desc(logLik), .by_group=T) %>% 
  slice(1) -> sums_species

# add absolute latitude
sums_species$abslatitude = abs(Site_table$lat[match(sums_species$site, Site_table$ID)])

# add species abundances
sums_species$abundance = NA
for (site in unique(sums_species$site)) {
  
  load(paste0(path_output, "meta_abundances/", site, "_abundances.Rdata"))
  
  # add abundances per species (N)
  sums_species$abundance[sums_species$site == site] = abundances$Nha[match(sums_species$sp[sums_species$site == site]
                                                                            , abundances$sp)]  
  
}

# plot latitudinal patterns in optimal decay_con
pdf(paste0(path_mortality, "/gridsearch_summary/optLL_con_species.pdf")
    , height = 4, width = 9)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), las = 1)

smoothScatter(sums_species$abslatitude, sums_species$decay_con
              , xlab = "latitude (°)"
              , ylab = expression(paste("optimal ", mu, " con"))
)

m = mgcv::gam(decay_con ~ s(abslatitude), data = sums_species)
x = seq(0, 55, 1)
pred = (mgcv::predict.gam(m, se.fit = T, newdata = list(abslatitude = x)))
lines(x, pred$fit, col = "red", lwd = 2)
polygon(c(x, rev(x)), c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit)), col = rgb(1, 0.1, 0.1, 0.3), border = NA)


smoothScatter(log(sums_species$abundance), sums_species$decay_con
              , xlab = "abundance (N/ha)"
              , ylab = expression(paste("optimal ", mu, " con"))
              , xaxt = "n"
              )
axis(1, at = log(10^(-1:4)), labels = 10^(-1:4), las = 1)

m = mgcv::gam(decay_con ~ s(log(abundance)), data = sums_species)
x = seq(min(sums_species$abundance), max(sums_species$abundance), 1)
pred = (mgcv::predict.gam(m, se.fit = T, newdata = list(abundance = x)))
lines(log(x), pred$fit, col = "red", lwd = 1.5)
polygon(c(log(x), rev(log(x)))
        , c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit))
        , col = rgb(1, 0.1, 0.1, 0.2), border = NA)

text(grconvertX(0.05, from = "ndc")
     , grconvertY(0.95, from = "ndc")
     , labels = "a"
     , cex = 3/2
     , xpd=NA)
text(grconvertX(0.5, from = "ndc")
     , grconvertY(0.95, from = "ndc")
     , labels = "b"
     , cex = 3/2
     , xpd=NA)

dev.off()




# Clean -------------------------------------------------------------------


rm(list = ls()[!(grepl("site", ls()) | grepl("path", ls()))])





