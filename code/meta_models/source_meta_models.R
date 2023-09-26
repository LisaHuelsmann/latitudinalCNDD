

# Fit meta-regressions and create outputs

library(gtsummary)
library(flextable)
library(officer)
library(dplyr)
library(tidyverse)
library(qgam)
library(metafor)
library(optimParallel)
library(purrr)
library(effects)
library(Hmisc)
library(sfsmisc)
library(sf)
library(spData)
library(spdep)
library(wesanderson)



# get path objects from paths
for (i in names(paths)) assign(i, paths[i])

# create folder for output
dir.create(paste0(path_meta, run))




# Load global mortality runs ----------------------------------------------


# load combined results from run
load(paste0(path_mortality, run, "/global_mortality.Rdata"))


# chose available results
terms = unique(AMEsamples_global$term)[startsWith(unique(AMEsamples_global$term), "con_")]
types = c("AME", "rAME")
changes = unique(AMEsamples_global$change)
changes = changes[changes != "iqr"]          # chose changes but exclude iqr (estimates not comparable!)
sites = sort(unique(AMEsamples_global$site))




# Settings ----------------------------------------------------------------



# suppress scientific notation
options(scipen = 9999)


# reference values for standardizing predictors in meta regression
ref_abund = 1
ref_lat = 11.75


# latitudes for plotting
latitude_names = c("Tropical", "Subtropical", "Temperate")
latitudes = c(11.75, 29.25, 45)
lat_breaks = seq(0, 52, 1)

# abundances for plotting
abundance = c(1, 10, 100)
# approximately 0.05, 0.5 and 0.95 quantiles



# Transformations ---------------------------------------------------------


# functions to transform response variables

# for AME, based on residuals
# trans_AME = function(x) sign(x)*(abs(x)^(1/3))
# backtrans_AME = function(x) sign(x)*((abs(x))^3)
trans_AME = function(x) x
backtrans_AME = function(x) x


# for rAME, based on common practice for RR
trans_rAME = function(x) log(x+1)
backtrans_rAME = function(x) exp(x)-1






# Uncertainty and significance --------------------------------------------


# point estimates from MLE
# significance and std.errors from samples (transformed)

for (type in types) {
  
  # get sample and aggregated object with MLE
  temp = get(paste0(type, "samples_global"))
  
  # get transformation
  trans = get(paste0("trans_", type))
  
  temp %>% 
    group_by(across(-c(estimate, MLE))) %>%     # group by all except estimate and MLE, results in grouping per site_sp x change with n = iter
    mutate(q2.5 = stats::quantile(estimate, 0.025),    #  quantiles for significance
           q97.5 = stats::quantile(estimate, 0.975),
           transSamples = trans(estimate),
           transMLE = trans(MLE)) %>% 
    summarise(std.error = diff(stats::quantile(transSamples, c(pnorm(-1), pnorm(1))))/2,
              estimate = unique(transMLE),
              significant = unique(q2.5 > 0 | q97.5 < 0)) %>% 
    as.data.frame() -> temp
  
  
  # store site as factor
  temp$site = as.factor(temp$site)
  
  
  assign(paste0(type, "sums_global"), temp)
  
}




# Add meta information -----------------------------------------------------


output_objects = gsub("_global", "", ls()[grepl("_global", ls())])

source("code/meta_models/source_add_meta_information.R", local = T)




# Demographic rates and tradeoffs -----------------------------------------


# calculate/transform species-level characteristics and demographic rates
sums_global %>% 
  transmute(site = site,
            sp = sp,
            abs_latitude = abs(latitude),
            temperature = MAT,
            precipitation = MAP,
            log_abundance = log(abundance),
            log1p_growth = log1p(incr_median),
            logit_survival = qlogis(1-mort_rate),
            log_size = log(dbh_q90)) -> M
M$log1p_growth[is.na(M$log1p_growth)] = 0

# calculate standardized demographic rates per site
M %>% 
  group_by(site) %>%
  mutate(stand_growth = as.vector(scale(log1p_growth)),
         stand_surv = as.vector(scale(logit_survival)),
         stand_size = as.vector(scale(log_size))) -> M_standardized


# PCA and demographic tradeoff axes via rotation
PCA = prcomp(M_standardized[, c("stand_growth", "stand_surv", "stand_size")], scale = TRUE)

r = 45 # chosen to mimic growth-mortality vs stature-recruitment tradeoff by Rüger et al.
coor = PCA$x[, 1:2]
coor_new <- Rotation(coor, r*pi/180)
M_standardized$tradeoff1 = coor_new[, 1]
M_standardized$tradeoff2 = coor_new[, 2]
PCs = PCA$rotation[, 1:2]
PCs_new <- Rotation(PCs, r*pi/180)

# plot PCA
pdf(paste0(path_meta, run, "/demographic_tradeoffs.pdf")
    , height = 5, width = 8)
par(mfrow = c(1, 2))
biplot(PCA, scale = F, cex = c(0.2, 1))
plot(M_standardized[, c("tradeoff1", "tradeoff2")]
     , pch = 16, cex = 0.5, col = "grey70", xlim = c(-8, 8), ylim = c(-8, 8)
     , xlab = "tradeoff1 (growth-mortality)"
     , ylab = "tradeoff2 (stature-recruitment)")
abline(h = 0, v = 0, col = "grey70")
arrows(x0 = c(0), y0 = 0, x1 = 4*PCs_new[, 1], y1 = 4*PCs_new[, 2])
dev.off()


# add standardized demographic rates and tradeoffs to (r)AME
for (type in types) {
  
  # get sample and aggregated object with MLE
  temp = get(paste0(type, "sums_global"))
  
  temp = cbind(temp, M_standardized[match(paste(temp$site, temp$sp), 
                                          paste(M_standardized$site, M_standardized$sp)),
                                    c("stand_growth", "stand_surv", "tradeoff1", "tradeoff2")])
  
  assign(paste0(type, "sums_global"), temp)
  
}









# Plot overview predictors and species characteristics --------------------


pdf(paste0(path_meta, run, "/data_overview.pdf")
    , height = 130/25.4
    , width = 136/25.4   # column-and-a-half
    , pointsize = 7
)

# plot abundance vs demography
par(mfrow = c(2, 2))
for (i in c("logit_survival", "log1p_growth")) {
  smoothScatter(M$log_abundance, M[, i]
                , xlab = "abundance (N/ha)", ylab = gsub("_", " ", i)
                , transformation = function(x) x^.5
                , xaxt = "n"
                , nrpoints = 0
                )
  axis(1, at = log(10^(-1:4)), labels = 10^(-1:4), las = 1)
  
  for (qu in c(0.25, 0.5, 0.75)) {
    form = formula(paste(i, "~ s(log_abundance, k = 5)"))
    m = qgam::qgam(form, data = M, qu = qu)
    
    pred = cbind(x = sort(M$log_abundance)
                 , as.data.frame(mgcv::predict.gam(m, se.fit = T))[order(M$log_abundance), ])
    lines(pred[, c("x", "fit")], col = "red", lwd = 1.5)
    polygon(c(pred$x, rev(pred$x))
            , c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit))
            , col = rgb(1, 0.1, 0.1, 0.2), border = NA)  
  }
}

text(grconvertX(0.02, from = "ndc")
     , grconvertY(0.98, from = "ndc")
     , labels = "a"
     , cex = 8/7
     , font = 2
     , xpd=NA)
text(grconvertX(0.5, from = "ndc")
     , grconvertY(0.98, from = "ndc")
     , labels = "b"
     , cex = 8/7
     , font = 2
     , xpd=NA)


# plot abundance vs latitude
smoothScatter(M$abs_latitude, M$log_abundance
              , ylab = "abundance (N/ha)", xlab = "latitude (°)"
              , transformation = function(x) x^.5
              , yaxt = "n"
              , nrpoints = 0
              )
axis(2, at = log(10^(-1:4)), labels = 10^(-1:4), las = 1)


for (qu in c(0.25, 0.5, 0.75)) {
  
  m = qgam::qgam(log_abundance ~ s(abs_latitude, k = 5), data = M, qu = qu)
  
  x = seq(min(M$abs_latitude), max(M$abs_latitude), 1)
  pred = mgcv::predict.gam(m, se.fit = T, newdata = list(abs_latitude = x))
  
  lines(x, pred$fit, col = "red", lwd = 1.5)
  polygon(c(x, rev(x))
          , c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit))
          , col = rgb(1, 0.1, 0.1, 0.2), border = NA)
}

text(grconvertX(0.02, from = "ndc")
     , grconvertY(0.45, from = "ndc")
     , labels = "c"
     , cex = 8/7
     , font = 2
     , xpd=NA)

dev.off()







# Figure 1 ----------------------------------------------------------------




# Subpanel order
if (length(sites) == 23) {
  
  site_order = rev(c("hkk", "fus", "mos", "lam", "kha", "pas"
                     , "sin", "edo", "len", "kor", "idc"
                     , "ama", "lpl", "bci", "luq", "ucsc", "wfdp", "scbi"
                     , "ldw", "wab", "serc", "wyw", "zof"))
  
} else { 
  
  site_order = sites
  
}

# Meta info sites
Site_table = readxl::read_xlsx(paste0(path_input, "global_metainfo/plot_sites_information.xlsx"), sheet = 1)

Site_table %>%
  dplyr::filter(ID %in% site_order) %>%
  mutate(cuts = cut(abs(lat), breaks = lat_breaks, right = F)) %>%
  arrange(match(ID, site_order)) -> Site_table
site_names = Site_table$site

# Color choice
cuts_levels = levels(Site_table$cuts)
cols = rev(wesanderson::wes_palette("Zissou1", length(cuts_levels), type = "continuous"))
cols_order = cols[match(Site_table$cuts, cuts_levels)]  # reorder to match site_order


# Convert to simple feature
sites_sf = Site_table %>%
  st_as_sf(coords = c("long", "lat"),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("+proj=wintri")


# Worldmap
dunia = spData::world %>%
  st_as_sf() %>%
  st_transform("+proj=wintri")

grat = st_graticule(lat = c(-90, -45, 0, 45, 90)
                    , lon = T) %>%
  st_transform("+proj=wintri")



pdf(paste0(path_meta, run, "/Fig1.pdf")
    , height = 247/25.4*0.7
    , width = 183/25.4)


site_models = list()
for (term in terms) {
  
  for (change in changes) {
    
    for (type in types) {
      
      # select output
      output = get(paste0(type, "sums_global"))
      sel = output$term == term & output$change == change
      dat = output[sel, ]
      
      # Layout
      par(oma = c(0.5, 0.6, 0, 0), fig = c(0, 1, 0, 1), new = F)
      plot.new()
      
      # axes labels
      mtext("abundance (N/ha)", side = 1, line = -0.5, outer = T, cex = 0.6)
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      mtext(paste("stabilizing CNDD", unit), side = 2, line = -0.1, outer = T, cex = 0.6)
      
      # circular coordinates
      multiCoords = circular_coor(order = site_order, inner = ifelse(length(site_order) == 23, 7, 5))
      multiCenter = cbind(rowMeans(multiCoords[, 1:2]), rowMeans(multiCoords[, 3:4]))
      
      # add map in center
      par(mar = c(0, 0, 0, 0), fig = c(0.34, 0.66, 0.35, 0.65), new = T)
      plot(st_geometry(dunia)
           , col = "grey90"
           , border = F
           , bg = NULL
           , xlim = c(st_bbox(dunia)[1]+5000000, st_bbox(dunia)[3]-1500000)
           , ylim = c(st_bbox(dunia)[2]+3000000, st_bbox(dunia)[4]))
      plot(st_geometry(grat)[2:4], add = T, col = "grey90", lwd = 2, lty = 2)
      
      # connect panels and site coordinates
      par(xpd = NA, new = T)
      for (i in 1:length(site_order)) {
        coor_temp = st_coordinates(sites_sf)[sites_sf$ID == site_order[i], ]
        lines(c(grconvertX(multiCenter[i, 1], from = "ndc"), coor_temp[1])
              , c(grconvertY(multiCenter[i, 2], from = "ndc"), coor_temp[2])
              , col = "grey80"
              , lwd = 0.8
        )
        
      }
      par(xpd = F, new = T)
      
      # add forest site coordinates
      plot(st_geometry(sites_sf), add = T, pch = 16, col = cols_order, cex = 0.6)
      
      # plot estimates versus abundance
      par(mar = c(1, 1.2, 0.8, 0.1))
      temp = plot_estimates_vsAbundance(dat, type = type, order = site_order, names = site_names
                                        , trans = get(paste0("trans_", type))
                                        , backtrans = get(paste0("backtrans_", type))
                                        , multiCoords = multiCoords
                                        , cols = cols_order
                                        , color.axes = "black"
                                        , markRare = T
                                        , returnMod = T
                                        )
      temp = temp[order(names(temp))]
      site_models[[paste(type, change)]] = as.data.frame(bind_rows(lapply(temp, function(x) try(broom::tidy(x))), .id = "site"))
      
    }
  }
}

dev.off()


save(site_models
     , file = paste0(path_meta, run, "/sitemodels.Rdata"))

rm(site_models)




# Fit site specific meta-regressions without moderators -----------------------

# Note that estimates are not back-transformed!


# Fit models and store in list
sitemean_list = list()
i = 0

for (term in terms) {
  
  for (change in changes) {
    
    for (type in types) {
      
      (i = i+1)
      
      # select output
      output = get(paste0(type, "sums_global"))
      sel = output$term == term & output$change == change
      dat_sel = output[sel, ]
      
      
      out = data.frame()
      
      # site-specific
      for (site in site_order) {
        
        # prepare escalc
        dat = escalc(measure = "GEN"
                     , yi = estimate
                     , sei = std.error
                     , data = dat_sel[dat_sel$site == site, ])
        
        
        # fit model
        mod = rma_fit(dat, moderator = ~ 1)
        
        # check for non-reliable results 
        check = rma_convergence(mod)
        mod = check$mod
        reliable = check$reliable
        
        
        # export estimates when reliable (all not back-transformed)
        out = rbind(out, data.frame(site = site
                                    , mean = ifelse(reliable, as.vector(mod$b), NA)
                                    , ci.lb = ifelse(reliable, mod$ci.lb, NA) 
                                    , ci.ub = ifelse(reliable,  mod$ci.ub, NA)
                                    , sigma = ifelse(reliable, sqrt(mod$tau2), NA)
                                    , latitude = unique(dat$latitude[dat$site == site])))
        
      }
      out$col = cols_order
      
      sitemean_list[[i]] = list("out" = out
                                , "type" = type
                                , "term" = term
                                , "change" = change)
    }
  }
}

# Return summary of site-specific estimates
sink(paste0(path_meta, run, "/site-specific.txt"))

for (i in 1:length(sitemean_list)) {
 
  change = sitemean_list[[i]]$change
  type = sitemean_list[[i]]$type
  out = sitemean_list[[i]]$out
  
  # run
  cat(change, type, "\n")  
  
  # sites with NA CNDD
  cat("sites with NA CNDD:\n")
  cat("n =", sum(is.na(out$mean)), " ")
  cat(out$site[which(is.na(out$mean))], "\n")  
  out = out[which(!is.na(out$mean)), ] # only non NA
  
  # sites with negative average CNDD (works also for transformed values log1p(0) = 0)
  cat("sites with negative average CNDD:\n")
  cat("n =", sum(out$mean < 0), " ")
  cat(out$site[which(out$mean < 0)], "\n")  
  out = out[which(out$mean > 0), ] # only positive CNDD
  
  # sites with CV < 0.4
  cat("sites with CV < 0.4:\n")
  out$CV = out$sigma/out$mean
  cat("n =", sum(out$CV < 0.4), " ")  
  cat(out$site[which(out$CV < 0.4)], "\n")
  
  cat("\n")
}
sink()






# Figure 4 ----------------------------------------------------------------


# Prediction 3

pdf(paste0(path_meta, run, "/Fig4.pdf")
    , height = 100/25.4
    , width = 183/25.4  # double column 
    , pointsize = 7
)


for (i in 1:length(sitemean_list)) {
  
  term = sitemean_list[[i]]$term
  change = sitemean_list[[i]]$change
  type = sitemean_list[[i]]$type
  out = sitemean_list[[i]]$out
  
  # calculate CV (transformed scale)
  out$CV = out$sigma/out$mean
  
  # take out values with negative average CNDD
  out = out[which(out$mean > 0), ]
  
  
  par(mfrow = c(1, 2)
      , mar = c(4, 5, 3, 1))
  
  # a) plot coefficient of variation against latitude
  plot(CV ~ abs(latitude)
       , out
       , las = 1
       , ylim = c(0, max(out$CV))
       , xlab = "absolute latitude (°)"
       , ylab = "coefficient of variation of CNDD"
       , bty = 'l'
       , type = 'n'
  )
  fit = lm(CV ~ abs(latitude), out)
  lats = seq(0, max(out$latitude))
  pred = predict(fit, newdata = data.frame(latitude = lats),  se.fit = T)
  lines(lats, pred$fit)
  polygon(c(lats, rev(lats))
          , c(pred$fit + 1.96*pred$se.fit, rev(pred$fit - 1.96*pred$se.fit))
          , col = add.alpha("black", 0.05)
          , border = F
  )
  abline(h = 0.4
         , col = "grey70"
         , lty = 3
         , lwd = 2
  )
  points(CV ~ abs(latitude), out
         , col = out$col
         , pch = 16
  )
  text(lats[1]
       , pred$fit[1] + 0.15
       , paste0("p=", round(summary(fit)$coef[2, 4], 2))
       , adj = 0
       , cex = 6/7
  )
  
  
  # b) plot standard deviation against mean (transformed scale)
  par(mar = c(4, 5, 3, 1))
  plot(sigma ~ mean
       , out
       , xlim = range(c(out$ci.lb, out$ci.ub))
       , ylim = c(0, max(out$sigma))
       , xlab = "mean CNDD (transformed scale)"
       , ylab = ""
       , bty = 'l'
       , las = 1
       , type = 'n'
  )
  mtext("standard deviation of CNDD (transformed scale)", 2, line = 4)
  
  # Lines with same CV
  for (cv in c(0.2, 0.4, 0.6)) {
    lines(seq(0, max(out$ci.ub), len = 10)
          , cv*seq(0, max(out$ci.ub), len = 10)
          , col = "grey70")
    text(max(out$ci.ub)
         , cv*max(out$ci.ub) + 0.01*max(out$sigma)
         , paste0("CV = ", cv)
         , col = "grey70"
         , adj = c(1, 0)
         , cex = 6/7)
  }
  
  # Means and sigmas of CNDD estimates
  segments(out$ci.lb, out$sigma, out$ci.ub, out$sigma, col = out$col)
  points(sigma ~ mean
         , out
         , col = out$col
         , pch = 16
  )
  
  # add label for panels
  text(grconvertX(c(0.02, 0.5), from = "ndc")
       , grconvertY(c(0.98, 0.98), from = "ndc")
       , labels = c("a", "b")
       , cex = 8/7
       , font = 2
       , xpd=NA)
  
}

dev.off()





# Autocorrelation summary -------------------------------------------------


if (any(!is.null(sums_global$test.spatial))) {
  
  collect_significance = function(x) {
    c(n = sum(x), perc = round(100*sum(x)/length(x), 2))
  }
  
  # p-values
  sums_global %>% 
    mutate(none = test.spatial < 0.05,
           holm = p.adjust(test.spatial) < 0.05) %>% 
    summarise(none = collect_significance(none),
              holm = collect_significance(holm)) %>% 
    t() -> test.spatial

  test.spatial = as.data.frame(test.spatial)
  names(test.spatial) = c("N", "%")
  test.spatial$adjustment = c("none", "holm")
  test.spatial = test.spatial[, c("adjustment", "N", "%")]
  
  
  sink(paste0(path_meta, run, "/autocorrelation.txt"))
  cat(names(test.spatial), "\n")
  cat(paste(test.spatial[1, ]), "\n")
  cat(paste(test.spatial[2, ]), "\n")
  sink()

}







# Fit global CNDD models ---------------------------------------------------------



# Fit models and store in list
mod_list = list()
i = 0

for (term in terms) {

  for (change in changes) {

    for (type in types) {

      (i = i+1)

      
      # select output
      output = get(paste0(type, "sums_global"))
      sel = output$term == term & output$change == change
      dat_sel = output[sel, ]


      # transformation of abundance and latitude
      transAbund = function(x, ref_abund) (log(x)-log(ref_abund))  # test for main effects are evaluated at abundance = ref_abund
      transLat = function(x, ref_lat)   (abs(x)-ref_lat)             #                                     and latitude = ref_lat
      dat_sel$tAbundance = transAbund(dat_sel$abundance, ref_abund)
      dat_sel$tLatitude = transLat(dat_sel$latitude, ref_lat)
      
      # prepare escalc
      dat = escalc(measure = "GEN"
                   , yi = estimate
                   , sei = std.error
                   , slab = paste(sp, site, sep = ", ")
                   , data = dat_sel)

      # fit abundance-mediated model
      mod = rma.mv(yi = yi
                   , V = vi
                   , mods = ~ tLatitude * tAbundance
                   , random = list( ~ 1 | site/sp)           # random intercept for species in site
                   , method = "REML"
                   , data = dat
                   , control=list(optimizer = "optimParallel", ncpus = ncpu_meta)
                   , sparse = T
                   # , verbose = T
      )


      # starting values (broad guesses to speed up fitting)
      start_sigma = mod$sigma2
      if (any(start_sigma <= 0)) start_sigma = rep(coef(mod)[1]/2, length(start_sigma))
      
      
      # fit mean species model
      mod1 = rma.mv(yi = yi
                    , V = vi
                    , mods = ~ tLatitude
                    , random = list(~ 1 | site/sp)
                    , method = "REML"
                    , data = dat
                    , control=list(optimizer = "optimParallel"
                                   , ncpus = ncpu_meta
                                   , sigma2.init = start_sigma)
                    , sparse = T
                    # , verbose = T
      )
    
      
      # reduced data models only for main run
      if (run == "main") {
        
        # Cooks distance
        # how much the whole regression model would change if ith case is removed
        D = cooks.distance(mod, reestimate = F, ncpus = ncpu_meta, parallel = "snow", progbar = F)
        exclude = D>0.005
        dat_reduced = dat[!(exclude | is.na(exclude)), ]
        sites_reduced = unique(dat_reduced$site)
        
        # fit abundance-mediated model with reduced data
        mod_x = rma.mv(yi = yi
                       , V = vi
                       , mods = ~ tLatitude * tAbundance
                       , random = list(~ 1 | site/sp)
                       , method = "REML"
                       , data = dat_reduced
                       , control=list(optimizer = "optimParallel"
                                      , ncpus = ncpu_meta
                                      , sigma2.init = start_sigma)
                       , sparse = T
                       # , verbose = T
        )
        
        # fit mean species model with reduced data
        mod1_x = rma.mv(yi = yi
                        , V = vi
                        , mods = ~ tLatitude
                        , random = list(~ 1 | site/sp)
                        , method = "REML"
                        , data = dat_reduced
                        , control=list(optimizer = "optimParallel"
                                       , ncpus = ncpu_meta
                                       , sigma2.init = start_sigma)
                        , sparse = T
                        # , verbose = T
        )
        
      } else {
        
        mod_x = mod1_x = NA
        sites_reduced = NA
        
      }
      
      
      
      # fit model with tradeoff axes
      mod2 = rma.mv(yi = yi
                    , V = vi
                    , mods = ~ tLatitude * tAbundance + tradeoff1 * tradeoff2
                    , random = list(~ 1 | site/sp)
                    , method = "REML"
                    , data = dat
                    , control=list(optimizer = "optimParallel"
                                   , ncpus = ncpu_meta
                                   , sigma2.init = start_sigma)
                    , sparse = T
                    # , verbose = T
      )


      # fit model with demographic rates
      mod3 = rma.mv(yi = yi
                    , V = vi
                    , mods = ~ tLatitude * tAbundance + (stand_growth + stand_surv)^2
                    , random = list(~ 1 | site/sp)
                    , method = "REML"
                    , data = dat
                    , control=list(optimizer = "optimParallel"
                                   , ncpus = ncpu_meta
                                   , sigma2.init = start_sigma)
                    , sparse = T
                    # , verbose = T
      )
     
      
      # create output
      out = list("mod0" = mod
                 , "mod1" = mod1
                 , "mod2" = mod2
                 , "mod3" = mod3
                 , "mod0x" = mod_x
                 , "mod1x" = mod1_x
                 , "data" = dat
                 , "type" = type
                 , "term" = term
                 , "change" = change
                 , "transAbund" = transAbund
                 , "transLat" = transLat
                 , "ref_abund" = ref_abund
                 , "ref_lat" = ref_lat
                 , "sites_reduced" = sites_reduced
                 )

      mod_list[[i]] = out

    }
  }
}

# save main models as list to use for comparisons
mod_list_red = lapply(mod_list, function(x) x[!names(x) %in% c("mod0x", "mod1x", "mod2", "mod3")])
save(list = c("mod_list_red", "ref_abund", "ref_lat", ls()[sapply(ls(), function(x) is.function(get(x)))])
     , file = paste0(path_meta, run, "/metamodels.Rdata"))


# wrap mod_list in res[[run]] 
# this allows that all plotting functions work for both a single and several runs
res = list()
res[[run]] = mod_list

rm(mod_list_red, mod_list)







# Plot joint meta-regression results ----------------------------------------------------------------


# a) Prediction 1 CNDD against latitude
# b) Prediction 2 CNDD against abundance at three latitudes: in color and at specified lats
# c) Prediction 2 CNDD against latitude for three abundances


# full dataset
pdf(paste0(path_meta, run, "/Meta_regressions.pdf")
    , height = 190/25.4
    , width = 136/25.4   # column-and-a-half
    , pointsize = 7
)

# loop through CNDD definitions
for (i in 1:length(res[[run]])) {

  # get results
  res_i = lapply(res, "[[", i)

  # back transformation
  backtrans = get(paste0("backtrans_", res_i[[1]]$type))
  
  # generate predictions
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1", pvalue = T)
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0", pvalue = T)
  preds_lat_abund = lapply(res_i, get_predictions_latitude, abundances = abundance, select = "mod0", pvalue = T)
  
  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_lat_abund, function(x) x[[1]]$ylim)),
                 backtrans(sitemean_list[[i]]$out$mean)*100), na.rm = T)

  # plotting
  layout(matrix(c(1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7), ncol = 3, byrow = T))
  par(las = 1
      , oma = c(3, 4, 0, 0)
      , mar = c(2, 2, 3, 1)
      , cex = 1)

  # only latitude model
  plot_latitude(preds_lat
                , labelsx = 1, labelsy = 1
                , ylim = ylim
                , col = "black")
  
  # add site-specific means
  out_trans = 100*backtrans(sitemean_list[[i]]$out[, c("mean", "sigma", "ci.lb", "ci.ub")])
  out_trans = cbind(out_trans, sitemean_list[[i]]$out[, c("site", "latitude", "col")])
  points(abs(out_trans$latitude)
         , out_trans$mean
         , pch = 16
         , col = add.alpha("black", 0.2))
  
  # latitude*abundance model latitude-moderated
  col = cols[match(cut(latitudes, breaks = lat_breaks, right = F), cuts_levels)]
  plot_abundance(preds_abund
                 , latitude_names = latitude_names
                 , labelsx = 2, labelsy = 1
                 , ylim = ylim
                 , col_lat = col)
  
  # latitude*abundance model abundance-moderated 
  plot_latitude(preds_lat_abund
                , labelsx = 2, labelsy = 1
                , ylim = ylim
                , col = "black")
  
  # add label for panels
  text(grconvertX(c(0.02, 0.02, 0.02), from = "ndc")
       , grconvertY(c(0.98, 0.48, 0.24), from = "ndc")
       , labels = c("a", "b", "c")
       , cex = 8/7
       , font = 2
       , xpd=NA)
  
}
dev.off()




# reduced dataset (only for main run)
if (run == "main") {

  pdf(paste0(path_meta, run, "/Meta_regressions_reduced_data.pdf")
      , height = 190/25.4
      , width = 136/25.4   # column-and-a-half
      , pointsize = 7
  )
  
  # loop through CNDD definitions
  for (i in 1:length(res[[run]])) {
    
    # get results
    res_i = lapply(res, "[[", i)
    
    # generate predictions
    preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1x", pvalue = T)
    preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0x", pvalue = T)
    preds_lat_abund = lapply(res_i, get_predictions_latitude, abundances = abundance, select = "mod0x", pvalue = T)
    
    # combine y limits
    ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                   unlist(lapply(preds_abund, function(x) x[[1]]$ylim)),
                   unlist(lapply(preds_lat_abund, function(x) x[[1]]$ylim)),
                   backtrans(sitemean_list[[i]]$out$mean)*100), na.rm = T)
    
    # plotting
    layout(matrix(c(1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7), ncol = 3, byrow = T))
    par(las = 1
        , oma = c(3, 4, 0, 0)
        , mar = c(2, 2, 3, 1)
        , cex = 1)
    
    # only latitude model
    plot_latitude(preds_lat
                  , labelsx = 1, labelsy = 1
                  , ylim = ylim
                  , col = "black")
    
    # latitude*abundance model latitude-moderated
    col = cols[match(cut(latitudes, breaks = lat_breaks, right = F), cuts_levels)]
    plot_abundance(preds_abund
                   , latitude_names = latitude_names
                   , labelsx = 2, labelsy = 1
                   , ylim = ylim
                   , col_lat = col)
    
    # latitude*abundance model abundance-moderated 
    plot_latitude(preds_lat_abund
                  , labelsx = 2, labelsy = 1
                  , ylim = ylim
                  , col = "black")
    
    # add label for panels
    text(grconvertX(c(0.02, 0.02, 0.02), from = "ndc")
         , grconvertY(c(0.98, 0.48, 0.24), from = "ndc")
         , labels = c("a", "b", "c")
         , cex = 8/7
         , font = 2
         , xpd=NA)
    
  }
  dev.off()
  
}





# Figure 2 ----------------------------------------------------------------

# Prediction 1

pdf(paste0(path_meta, run, "/Fig2.pdf")
    , height = 110/25.4
    , width = 136/25.4   # column-and-a-half
    , pointsize = 7
    )

# loop through CNDD definitions
for (i in 1:length(res[[run]])) {
  
  # get results
  res_i = lapply(res, "[[", i)
  
  # back transformation
  backtrans = get(paste0("backtrans_", res_i[[1]]$type))
  
  # generate predictions
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1", pvalue = T)
  
  # other predictions for y limits
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0", pvalue = F)
  preds_lat_abund = lapply(res_i, get_predictions_latitude, abundances = abundance, select = "mod0", pvalue = F)
  
  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_lat_abund, function(x) x[[1]]$ylim)),
                 backtrans(sitemean_list[[i]]$out$mean)*100), na.rm = T)
  
  # plotting
  par(las = 1, mar = c(4, 5, 1, 1))
  
  # only latitude model
  plot_latitude(preds_lat
                , labelsx = 1, labelsy = 1
                , ylim = ylim
                , col = "black"
                , plocation = "right"
  )
  
  # add site-specific means
  out_trans = 100*backtrans(sitemean_list[[i]]$out[, c("mean", "sigma", "ci.lb", "ci.ub")])
  out_trans = cbind(out_trans, sitemean_list[[i]]$out[, c("site", "latitude", "col")])
  
  # label sites
  out_trans$site_names = paste(" ", Site_table$site[match(out_trans$site, Site_table$ID)], " ")
  addTextLabels(abs(out_trans$latitude)
                , out_trans$mean
                , gsub('(.{1,30})(\\s|$)', '\\1\n', out_trans$site_names)
                , cex.label = 5/7
                , col.label = "grey40"
                , col.line = "grey70"
                , lwd = 0.8
                , cex.pt = 0.9
  )
  
  # dots
  points(abs(out_trans$latitude)
         , out_trans$mean
         , pch = 16
         , col = "grey40"
  )
  
  
}
dev.off()






# Figure 3 ----------------------------------------------------------------


# Prediction 2

# a) CNDD against abundance at three latitudes: in color and at specified lats
# b) CNDD against latitude for three abundances


pdf(paste0(path_meta, run, "/Fig3.pdf")
    , height = 130/25.4
    , width = 183/25.4  # double column 
    , pointsize = 7
)

# loop through CNDD definitions
for (i in 1:length(res[[run]])) {
  
  # get results
  res_i = lapply(res, "[[", i)
  
  # back transformation
  backtrans = get(paste0("backtrans_", res_i[[1]]$type))
  
  # generate predictions
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0", pvalue = T)
  preds_lat_abund = lapply(res_i, get_predictions_latitude, abundances = abundance, select = "mod0", pvalue = T)
  
  # other predictions for y limits
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1", pvalue = F)
  
  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_lat_abund, function(x) x[[1]]$ylim)),
                 backtrans(sitemean_list[[i]]$out$mean)*100), na.rm = T)
  
  # plotting
  layout(matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = T))
  par(las = 1
      , oma = c(3, 4, 0, 0)
      , mar = c(2, 2, 3, 1)
      , cex = 1)
  
  # latitude*abundance model latitude-moderated
  col = cols[match(cut(latitudes, breaks = lat_breaks, right = F), cuts_levels)]
  plot_abundance(preds_abund
                 , latitude_names = latitude_names
                 , labelsx = 2, labelsy = 1
                 , ylim = ylim
                 , col_lat = col)
  
  # latitude*abundance model abundance-moderated 
  plot_latitude(preds_lat_abund
                , labelsx = 2, labelsy = 1
                , ylim = ylim
                , col = "black")
  
  # add label for panels
  text(grconvertX(c(0.02, 0.02), from = "ndc")
       , grconvertY(c(0.98, 0.48), from = "ndc")
       , labels = c("a", "b")
       , cex = 8/7
       , font = 2
       , xpd = NA)
  
}
dev.off()



# Plot CNDD model diagnostics ---------------------------------------------


# https://stats.stackexchange.com/questions/155693/metafor-package-bias-and-sensitivity-diagnostics

pdf(paste0(path_meta, run, "/model_checks.pdf"))

lapply(res[[run]], function(x) {

  par(mfrow = c(2, 2), las = 1, mar = c(4, 2, 3, 1))

  # get residuals
  rs = residuals(x$mod0, type = "rstandard") # rstudent is very slow because of refitting the model
  d = as.data.frame(x$mod0$X)
  d$res = rs

  # distribution check
  stats::qqnorm(rs, col = rgb(0, 0, 0, 0.4), cex = 0.3)
  stats::qqline(rs)

  # against predicted
  d$fitted = fitted(x$mod0)
  plot(d$fitted, d$res)
  xlim = range(d$fitted)
  x_line = seq(xlim[1], xlim[2], length.out = 200)
  xs = data.frame(fitted = x_line)
  q = qgam::mqgam(res ~ s(fitted)
                  , qu = c(0.25, 0.5, 0.75), data = d)
  for(iq in c(0.25, 0.5, 0.75)){
    pred <- qdo(q, iq, predict, newdata = xs)
    lines(xs$fitted, pred, col = 2)
  }

  # against latitude
  plot(d$tLatitude, d$res)
  xlim = c(0, max(x$data$latitude))
  x_line = seq(xlim[1], xlim[2], length.out = 200)
  xs = data.frame(tLatitude = x$transLat(x_line, ref_lat))
  q = qgam::mqgam(res ~ s(tLatitude)
                  , qu = c(0.25, 0.5, 0.75), data = d)
  for(iq in c(0.25, 0.5, 0.75)){
    pred <- qdo(q, iq, predict, newdata = xs)
    lines(xs$tLatitude, pred, col = 2)
  }

  # against abundance
  plot(d$tAbundance, d$res)
  xlim = c(min(x$data$abundance), max(x$data$abundance))
  x_line = seq(xlim[1], xlim[2], length.out = 200)
  xs = data.frame(tAbundance = x$transAbund(x_line, ref_abund))
  q = qgam::mqgam(res ~ s(tAbundance)
                  , qu = c(0.25, 0.5, 0.75), data = d)
  for(iq in c(0.25, 0.5, 0.75)){
    pred <- qdo(q, iq, predict, newdata = xs)
    lines(xs$tAbundance, pred, col = 2)
  }

})

dev.off()





# Write CNDD model summaries ----------------------------------------------


if (run == "main") {
  
  # reduced dataset models only for main run
  settings = data.frame(names = c("mod0", "mod1", "mod2", "mod3", "mod0x", "mod1x")
                        , labels = c("main", "noabund", "tradeoffs", "demography", "main_reduced", "noabund_reduced"))
  
} else {
  settings = data.frame(names = c("mod0", "mod1", "mod2", "mod3")
                        , labels = c("main", "noabund", "tradeoffs", "demography"))
}


for (i in settings$names) {

  
  # write function
  write_word_table <- function(var, doc){
    doc %>%
      body_add_flextable(var) %>%
      body_add_par("") %>%
      body_add_par("") }

  # list of tables
  my_list <- lapply(res[[run]], function(x) {
    
    # digits depending on effect size
    temp = format(coef(x[[i]]), signdigits = 2, digits = 2, trim = T)
    digits = decimalplaces(temp)
    
    # collect model summary
    tbl_regression(x[[i]]
                   , estimate_fun = ~style_sigfig(., digits = digits)
                   , pvalue_fun = function(x) case_when(is.na(x) ~ NA_character_,
                                                        x < 0.0001 ~ format(x, digits = 2, scientific = TRUE),
                                                        x >= 0.0001 ~ formatC(signif(x, digits = 2), digits = 2))
                   ) %>%
      bold_p() %>%
      as_flex_table() %>%
      set_caption(caption = paste("CNDD assessed as", x$type, "in", x$term, "calculated at", x$change,
                                  "densities in %.", "Fitted with", class(x$mod0)[1])) %>%
      font(fontname = "Times New Roman", part = "all") %>%
      
      # standard deviations for site and species
      add_footer_lines(values = paste("sigma", x[[i]]$s.names, "=", 
                                      format(sqrt(x[[i]]$sigma2), nsmall = 2, digits = 2))) %>% 
      
      # pseudo R2 for including additional life history strategies
      # proportional reduction in the variance components as a sort of pseudo R-squared value
      # https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model-in-metafor
      add_footer_lines(values = ifelse(i %in% c("mod2", "mod3")
                                       , paste0("Pseudo R2 = ", 
                                                round((sum(x$mod0$sigma2) - sum(x[[i]]$sigma2)) / sum(x$mod0$sigma2), 3))
                                       , NA)) %>% 
      
      # add n
      add_footer_lines(values = paste("n", "=", sum(x[[i]]$not.na))) %>%
      
      # add n sites
      add_footer_lines(values = ifelse(i %in% c("mod0x", "mod1x")
                                       , paste0("N sites = ", length(x$sites_reduced))
                                       , paste0("N sites = ", 23))) %>% 
      
      theme_booktabs()
  })
  

  # use walk (the invisible function of map) to include all tables in one doc
  my_doc <- read_docx()
  walk(my_list, write_word_table, my_doc)
  print(my_doc, target = paste0(path_meta, run, "/model_summaries_"
                                , settings$labels[which(settings$names == i)], ".docx")) %>% invisible()

}



# IQR ---------------------------------------------------------------------


# collect means and CI
iqr_res = list()

# only when iqr result is available
if (sum(AMEsums_global$change == "iqr") > 0) {
  
  
  pdf(paste0(path_meta, run, "/iqr.pdf")
      , height = 100/25.4
      , width = 89/25.4  # single column 
      , pointsize = 7
  )
  par(mar = c(4, 5, 1, 1)
      , las = 1)
  
  change = "iqr"
  
  for (term in terms) {
    
    for (type in types) {
      
      # select output
      output = get(paste0(type, "sums_global"))
      sel = output$term == term & output$change == change
      dat_sel = output[sel, ]
      
      # backtransformation
      backtrans = get(paste0("backtrans_", type))
      
      # plotting settings
      xlim = c(-1, 1)*max(abs(backtrans(dat_sel$estimate))*100)/20
      if (type == "rAME") xlim[1] = -100
      if (type == "rAME") ylim = c(0, 500) else ylim = c(0, 800)
      breaks = ifelse(type == "rAME", 1000, 2000)
      outside = sum(backtrans(dat_sel$estimate)*100 < xlim[1] | backtrans(dat_sel$estimate)*100 > xlim[2]) / nrow(dat_sel)
      
      
      # histogram
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      hist(backtrans(dat_sel$estimate)*100
           , xlab = paste("stabilizing CNDD", unit)
           , xlim = xlim
           , main = ""
           , breaks = breaks
           , border = "grey70"
           , ylim = ylim)
      abline(v = 0)
      
      # add quantile range
      qs = quantile(backtrans(dat_sel$estimate)*100, probs = c(0.25, 0.75))
      arrows(qs[1], ylim[2]*0.8, qs[2], ylim[2]*0.8, angle = 90, code = 3, length = 0.04)
      interpretation = ifelse(type == "AME"
                              , "change in mortality per year"
                              , "relative change in mortality")
      text(qs[2]*2, ylim[2]*0.8, paste0("50 % of CNDD estimates between\n"
                                        , round(qs[1], 2), " and ", round(qs[2], 2)
                                        , " %\n", interpretation)
           , cex = 6/7, adj = 0)
      
      
      
      # add mean and CI from meta-regression
      dat = escalc(measure = "GEN"
                   , yi = estimate
                   , sei = std.error
                   , slab = paste(sp, site, sep = ", ")
                   , data = dat_sel)
      mod = rma.mv(yi = yi
                   , V = vi
                   , random = list( ~ 1 | site/sp)           # random intercept for species in site
                   , method = "REML"
                   , data = dat
                   , control=list(optimizer = "optimParallel", ncpus = ncpu_meta)
                   , sparse = T
                   # , verbose = T
      )

      CI = predict(mod, transf = backtrans)
      arrows(CI$ci.lb*100, ylim[2]*0.6, CI$ci.ub*100, ylim[2]*0.6, angle = 90, code = 3, length = 0.04)
      points(backtrans(coef(mod))*100, ylim[2]*0.6
             , pch = 18, cex = 1.4, col = "darkred")
      text(qs[2]*2, ylim[2]*0.6, paste0("mean CNDD from meta-regression \nat "
                                        , round(backtrans(coef(mod))*100, 2), " % ", interpretation)
           , cex = 6/7
           , adj = 0)
      iqr_res[[paste(term, type, sep = "_")]] = c(mean =  backtrans(coef(mod))*100
                                                  , ci.lb = CI$ci.lb*100
                                                  , ci.ub = CI$ci.ub*100)
      
      # number of estimates outside
      text(xlim[1], ylim[1]+0.92*ylim[2], paste0(round(100*outside, 0), "% of estimates outside")
           , cex = 5/7
           , adj = c(0, 0)
           , col = "grey70")
      
    }
  }
  
  dev.off()

  # Return mean and CI of iqr results
  sink(paste0(path_meta, run, "/iqr.txt"))
  
  for (i in 1:length(iqr_res)) {
    
    cat(names(iqr_res[i]), "\n")
    cat(names(iqr_res[[i]]), "\n")
    cat(iqr_res[[i]], "\n", "\n")
    c
    
  }
  sink()
  
    
}

