

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
sites = sort(unique(AMEsamples_global$site))




# Settings ----------------------------------------------------------------


# reference values for standardizing predictors in meta regression
ref_abund = 1
ref_lat = 11.75


# latitudes for plotting
latitude_names = c("Tropical", "Subtropical", "Temperate")
latitudes = c(11.75, 29.25, 45)
lat_breaks = seq(0, 52, 1)



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
            log_mortality = log(mort_rate),
            log_size = log(dbh_q90)) -> M

# calculate standardized demographic rates per site
M %>% 
  group_by(site) %>%
  mutate(stand_growth = as.vector(scale(log1p_growth)),
         stand_surv = as.vector(scale(qlogis(1-exp(log_mortality)))),
         stand_size = as.vector(scale(log_size))) -> M_standardized


# PCA and demographic tradeoff axes via rotation
PCA = prcomp(M_standardized[, c("stand_growth", "stand_surv", "stand_size")], scale = TRUE)

r = 45 # chosen to mimic growth-mortality vs stature-recruitment tradeoff by R??ger et al.
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
     , pch = 16, cex = 0.5, col = "grey", xlim = c(-8, 8), ylim = c(-8, 8)
     , xlab = "tradeoff1 (growth-mortality)"
     , ylab = "tradeoff2 (stature-recruitment)")
abline(h = 0, v = 0, col = "grey")
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
    , height = 8, width = 8)

# plot abundance vs demography
par(mfrow = c(2, 2))
for (i in c("log_mortality", "log1p_growth")) {
  smoothScatter(M$log_abundance, M[, i]
                , xlab = "abundance (N/ha)", ylab = gsub("_", " ", i)
                , transformation = function(x) x^.5
                , xaxt = "n"
                )
  axis(1, at = log(10^(-1:4)), labels = 10^(-1:4), las = 1)
  m = mgcv::gam(M[, i] ~ s(M$log_abundance))
  pred = cbind(x = sort(M$log_abundance)
               , as.data.frame(mgcv::predict.gam(m, se.fit = T))[order(M$log_abundance), ])
  lines(pred[, c("x", "fit")], col = "red", lwd = 2)
  polygon(c(pred$x, rev(pred$x)), c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit)), col = rgb(1, 0.1, 0.1, 0.3), border = NA)
}
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


# plot abundance vs latitude
smoothScatter(M$abs_latitude, M$log_abundance
              , ylab = "abundance (N/ha)", xlab = "latitude (??)"
              , transformation = function(x) x^.5
              , yaxt = "n"
              )
axis(2, at = log(10^(-1:4)), labels = 10^(-1:4), las = 1)
m = try(mgcv::gam(M$log_abundance ~ s(M$abs_latitude, k = 5)), silent = T)
text(grconvertX(0.05, from = "ndc")
     , grconvertY(0.45, from = "ndc")
     , labels = "c"
     , cex = 3/2
     , xpd=NA)

dev.off()







# Figure 1 ----------------------------------------------------------------




# Subpanel order
if (length(sites) == 23) {
  
  site_order = rev(c("hkk", "fus", "mos", "lam", "kha", "pas"
                     , "sin", "edo", "len", "idc", "kor"
                     , "ama", "lpl", "bci", "luq", "ucsc", "wfdp", "ldw"
                     , "wab", "scbi", "serc", "wyw", "zof"))
  
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
    , height = 165/25.4*0.9
    , width = 183/25.4)


for (term in terms) {
  
  for (change in changes) {
    
    for (type in types) {
      
      # select output
      output = get(paste0(type, "sums_global"))
      sel = output$term == term & output$change == change
      dat = output[sel, ]
      
      # Layout
      par(oma = c(0.5, 0.5, 0, 0), fig = c(0, 1, 0, 1), new = F)
      plot.new()
      
      # axes labels
      mtext("abundance (N/ha)", side = 1, line = -0.5, outer = T, cex = 0.6)
      mtext("stabilizing CNDD (%)", side = 2, line = -0.5, outer = T, cex = 0.6)
      
      # circular coordinates
      multiCoords = circular_coor(order = site_order, inner = 8)
      multiCenter = cbind(rowMeans(multiCoords[, 1:2]), rowMeans(multiCoords[, 3:4]))
      
      # add map in center
      par(mar = c(0, 0, 0, 0), fig = c(0.34, 0.66, 0.35, 0.65), new = T)
      plot(st_geometry(dunia)
           , col = "grey94"
           , border = F
           , bg = NULL
           , xlim = c(st_bbox(dunia)[1]+5000000, st_bbox(dunia)[3]-1500000)
           , ylim = c(st_bbox(dunia)[2]+3000000, st_bbox(dunia)[4]))
      plot(st_geometry(grat)[2:4], add = T, col = "grey94", lwd = 2, lty = 2)
      
      # connect panels and site coordinates
      par(xpd = NA, new = T)
      for (i in 1:length(site_order)) {
        coor_temp = st_coordinates(sites_sf)[sites_sf$ID == site_order[i], ]
        lines(c(grconvertX(multiCenter[i, 1], from = "ndc"), coor_temp[1])
              , c(grconvertY(multiCenter[i, 2], from = "ndc"), coor_temp[2])
              , col = "lightgrey" # cols_order[i]
              , lwd = 0.8
        )
        
      }
      par(xpd = F, new = T)
      
      # add forest site coordinates
      plot(st_geometry(sites_sf), add = T, pch = 16, col = cols_order, cex = 0.6)
      # text(st_coordinates(sites_sf), sites_sf$Abbreviation)
      
      # plot estimates versus abundance
      par(mar = c(1, 1.2, 1, 0.1))
      temp = plot_estimates_vsAbundance(dat, type = type, order = site_order, names = site_names
                                        , trans = get(paste0("trans_", type))
                                        , backtrans = get(paste0("backtrans_", type))
                                        , multiCoords = multiCoords
                                        , cols = cols_order, color.axes = "black" # "grey30"
                                        , markRare = T)
      
    }
  }
}

dev.off()








# Fit site specific meta-regressions without moderators -----------------------

# Note that estimates are untransformed!


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
        
        
        # export estimates when reliable (all untransformed)
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






# Figure 3 ----------------------------------------------------------------




pdf(paste0(path_meta, run, "/Fig3.pdf")
    , height = 5.5, width = 11)

for (i in 1:length(sitemean_list)) {
  
  term = sitemean_list[[i]]$term
  change = sitemean_list[[i]]$change
  type = sitemean_list[[i]]$type
  out = sitemean_list[[i]]$out
  
  # calculate CV (transformed scale)
  out$CV = out$sigma/out$mean
  
  # take out values with negative average CNDD
  out = out[which(out$mean > 0), ]
  
  
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  
  # a) plot coefficient of variation against latitude
  plot(CV ~ abs(latitude), out
       , col = out$col, pch = 16, las = 1
       , ylim = c(0, max(out$CV))
       , xlab = "absolute latitude (??)"
       , ylab = "coefficient of variation of CNDD"
       , bty='l'
  )
  fit = lm(CV ~ abs(latitude), out)
  lats = seq(0, max(out$latitude))
  pred = predict(fit, newdata = data.frame(latitude = lats),  se.fit = T)
  lines(lats, pred$fit)
  polygon(c(lats, rev(lats))
          , c(pred$fit + 1.96*pred$se.fit, rev(pred$fit - 1.96*pred$se.fit))
          , col = add.alpha("black", 0.05)
          , border = F)
  abline(h = 0.4, col = "grey", lty = 2)
  text(lats[1], pred$fit[1] + 0.1
       , paste0("p=", round(summary(fit)$coef[2, 4], 2))
       , adj = 0)
  
  
  # b) plot standard deviation against mean (transformed scale)
  par(mar = c(5, 6, 3, 1))
  plot(sigma ~ mean, out
       , col = out$col, pch = 16
       , xlim = range(c(out$ci.lb, out$ci.ub))
       , ylim = c(0, max(out$sigma))
       , xlab = "mean CNDD (transformed scale)"
       , ylab = ""
       , bty='l'
       , las = 1)
  mtext("standard deviation of CNDD (transformed scale)", 2, line = 4)
  segments(out$ci.lb, out$sigma, out$ci.ub, out$sigma, col = out$col)
  for (cv in c(0.2, 0.4, 0.6)) {
    lines(seq(0, max(out$ci.ub), len = 10), cv*seq(0, max(out$ci.ub), len = 10), col = "grey")
    text(0.9*max(out$ci.ub), 1.1*cv*max(out$sigma), paste0("CV = ", cv), col = "grey")
  }
  
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







# Fit CNDD models ---------------------------------------------------------



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
      
      # Cooks distance
      # how much the whole regression model would change if ith case is removed
      D = cooks.distance(mod, reestimate = F, ncpus = ncpu_meta, parallel = "snow", progbar = F)
      exclude = D>0.005
      dat_reduced = dat[!(exclude | is.na(exclude)), ]


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
                 )

      mod_list[[i]] = out

    }
  }
}

# save reduced model list to use for comparisons
mod_list_red = lapply(mod_list, function(x) x[!names(x) %in% c("mod0x", "mod1x", "mod2", "mod3")])
save(list = c("mod_list_red", "ref_abund", "ref_lat", ls()[sapply(ls(), function(x) is.function(get(x)))])
     , file = paste0(path_meta, run, "/metamodels.Rdata"))


# wrap mod_list in res[[run]] 
# this allows that all plotting functions work for both a single and several runs
res = list()
res[[run]] = mod_list

rm(mod_list_red, mod_list)







# Figure 2 ----------------------------------------------------------------


# 1) reduced model (without abundance)
# 2-4) CNDD against abundance at three latitudes: in color and at specified lats

# full dataset
pdf(paste0(path_meta, run, "/Fig2.pdf")
    , height = 7)

# loop through CNDD definitions
for (i in 1:length(res[[run]])) {

  # get results
  res_i = lapply(res, "[[", i)

  # back transformation
  backtrans = get(paste0("backtrans_", res_i[[1]]$type))
  
  # generate predictions
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1", pvalue = T)
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0", pvalue = T)

  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim)),
                 backtrans(sitemean_list[[i]]$out$mean)*100), na.rm = T)

  # plotting
  layout(matrix(c(1, 1, 1, 2, 3, 4), ncol = 3, byrow = T))
  par(las = 1, oma = c(3, 4, 0, 0), mar = c(2, 2, 3, 1))

  # only latitude model
  plot_latitude(preds_lat
                , labelsx = 1, labelsy = 1
                , ylim = ylim
                , col = "black"
                , panel = "a")
  
  # add site-specific means
  out_trans = 100*backtrans(sitemean_list[[i]]$out[, c("mean", "sigma", "ci.lb", "ci.ub")])
  out_trans = cbind(out_trans, sitemean_list[[i]]$out[, c("site", "latitude", "col")])
  points(abs(out_trans$latitude)
         , out_trans$mean
         , pch = 16
         , col = add.alpha("black", 0.2))
  
  # latitude*abundance model
  col = cols[match(cut(latitudes, breaks = lat_breaks, right = F), cuts_levels)]
  plot_abundance(preds_abund
                 , latitude_names = latitude_names
                 , labelsx = 2, labelsy = 1
                 , ylim = ylim
                 , col_lat = col
                 , panel = "b")

}
dev.off()


# reduced dataset
pdf(paste0(path_meta, run, "/Fig2_reduced_data.pdf")
    , height = 7)

# loop through CNDD definitions
for (i in 1:length(res[[run]])) {

  # get results
  res_i = lapply(res, "[[", i)

  # generate predictions
  preds_lat = lapply(res_i, get_predictions_latitude, abundances = 1,  select = "mod1x", pvalue = T)
  preds_abund = lapply(res_i, get_predictions_abundance, latitudes = latitudes, select = "mod0x", pvalue = T)

  # combine y limits
  ylim = range(c(unlist(lapply(preds_lat, function(x) x[[1]]$ylim)),
                 unlist(lapply(preds_abund, function(x) x[[1]]$ylim))))

  # plotting
  layout(matrix(c(1, 1, 1, 2, 3, 4), ncol = 3, byrow = T))
  par(las = 1, oma = c(3, 4, 0, 0), mar = c(2, 2, 3, 1))

  # only latitude model
  plot_latitude(preds_lat
                , labelsx = 1, labelsy = 1
                , ylim = ylim
                , col = "black"
                , panel = "a")

  # latitude*abundance model
  col = cols[match(cut(latitudes, breaks = lat_breaks, right = F), cuts_levels)]
  plot_abundance(preds_abund
                 , latitude_names = latitude_names
                 , labelsx = 2, labelsy = 1
                 , ylim = ylim
                 , col_lat = col
                 , panel = "b")

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



settings = data.frame(names = c("mod0", "mod1", "mod2", "mod3", "mod0x", "mod1x")
                      , labels = c("main", "noabund", "tradeoffs", "demography", "main_reduced", "noabund_reduced"))

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
    digits = ceiling(log10(1/min(abs(1*coef(x[[i]])))))
    
    # collect model summary
    tbl_regression(x[[i]]
                   , estimate_fun = ~style_sigfig(., digits = digits)
                   , pvalue_fun = function(x) scales::pvalue(x, accuracy = 0.01, add_p = F)) %>%
      bold_p() %>%
      # bold_labels() %>%
      as_flex_table() %>%
      set_caption(caption = paste("CNDD assessed as", x$type, "in", x$term, "calculated at", x$change,
                                  "densities in %.", "Fitted with", class(x$mod0)[1])) %>%
      font(fontname = "Times New Roman", part = "all") %>%
      
      # standard deviations for site and species
      add_footer_lines(values = paste("sigma", x[[i]]$s.names, "=", round(sqrt(x[[i]]$sigma2), digits))) %>% 
      
      # pseudo R2 for including additional life history strategies
      # proportional reduction in the variance components as a sort of pseudo R-squared value
      # https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model-in-metafor
      add_footer_lines(values = ifelse(i %in% c("mod2", "mod3")
                                       , paste0("Pseudo R2 = ", 
                                                round((sum(x$mod0$sigma2) - sum(x[[i]]$sigma2)) / sum(x$mod0$sigma2), 3))
                                       , NA)) %>% 
      
      # add n
      add_footer_lines(values = paste("n", "=", sum(x[[i]]$not.na))) %>%
      
      
      theme_booktabs()
  })
  
  my_doc <- read_docx()

  # use walk (the invisible function of map) to include all tables in one doc
  walk(my_list, write_word_table, my_doc)
  print(my_doc, target = paste0(path_meta, run, "/model_summaries_"
                                , settings$labels[which(settings$names == i)], ".docx")) %>% invisible()

}



