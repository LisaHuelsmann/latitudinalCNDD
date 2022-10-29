

# Mortality analyses
# Species-specific models: main





# Packages and data -------------------------------------------------------

library(dplyr)
library(tidyr)
library(broom)
library(mgcv)
library(mgcViz)
library(DHARMa)

# import helper functions
source("code/mortality_models/functions_marginal_effects_gam.R")

# get path objects from paths
for (i in names(paths)) assign(i, paths[i])





# Settings ----------------------------------------------------------------


# chose decay type
decay_type = "exp"

dec_fun <- function(sigma, distance, type) {  # define decay function explicitly

  if (type == "exp") {
    return( exp(-(1/sigma)*distance) )  
  }
  
  if (type == "expn") {
    return( exp(-(1/(sigma^2))*(distance^2)) )  
  }
}


# chose predictors for local density
decay_con = 3
decay_tot = 17
predictors = c(con_BA = paste0("ba_con_", decay_type, "_", decay_con),
               all_BA = paste0("ba_all_", decay_type, "_", decay_tot))



# settings for rare species handling
x = 20                                  # fixed cutoff for ndead and nsurv
nval = 4                                # minimum number of values for conspecific density

dbh_neighbor = 20                       # dbh of additional neighbor in mm
dist_neighbor = 1                       # distance to additional neighbor in m

minrange =                                        # minimum range for conspecific density:  
  1 *                                             # one time what is needed for marginal effects 
  pi*((dbh_neighbor/1000)/2)^2 *                  # basal area (m2) of one more neighbor with dbh_neighbor
  dec_fun(decay_con, dist_neighbor, decay_type)   # at a distance of dist_neighbor



# change in conspecific density for AME calculations
additive = pi*((dbh_neighbor/1000)/2)^2 *         # basal area (m2) of one more neighbor with dbh_neighbor
  dec_fun(decay_con, dist_neighbor, decay_type)   # at a distance of dist_neighbor  





# Load data ---------------------------------------------------------------


load(paste0(path_output, "data_3a_mortality/", site, "_tree_3a_mortality_", decay_type, ".Rdata"))
load(paste0(path_output, "meta_stature/", site, "_stature.Rdata"))


# rename data object
dat_mort = tree_mort
rm(tree_mort)


# generate predictor columns
for (i in 1:length(predictors)) {
  dat_mort[, names(predictors)[i]] = dat_mort[, predictors[i]]
}





# Data selection ----------------------------------------------------------


# remove NA in sp
dat_mort = dat_mort[!is.na(dat_mort$sp), ]

# remove ferns and palms
load(paste0(path_input, "data_species/", site, "_species.Rdata"))
dat_mort = dat_mort[dat_mort$sp %in% species$sp[species$fern.palm == "FALSE"], ]

# remove NA and 0 in interval
dat_mort = dat_mort[!is.na(dat_mort$interval) & !(dat_mort$interval == 0), ]

# limit to small trees
dat_mort = dat_mort[dat_mort$dbh < 100, ]

# save census as factor
dat_mort$census = as.factor(dat_mort$census)
table(dat_mort$census)




# Handling rare species ---------------------------------------------------



# explore available data per species
dat_mort %>% 
  group_by(sp) %>% 
  summarise(ndead = sum(surv_next==0),
            nsurv = sum(surv_next==1),
            range_con_BA = max(con_BA) - min(con_BA),
            max_con_BA = max(con_BA),
            unique_con_BA = length(unique(con_BA)),
            unique_all_BA = length(unique(all_BA)),
            unique_dbh = length(unique(dbh))
  ) %>% 
  
  # add issues
  mutate(issue_nval = unique_con_BA < nval,                        # less than nval unique values in conspecific densities
         issue_range = range_con_BA < minrange,                    # range should at least be equal to minrange
         issue_cutoff = ndead < x | nsurv < x,                     # ndead or nsurv below cutoff
         trymodel = !(issue_nval | issue_range | issue_cutoff),    # should a species specific model be tried at all?
         rare = !trymodel                                          # preliminary assignment of rare species
  ) -> nsp


# explore consequence of fixed cutoff of at least x observations of survival and death
plot(log1p(ndead) ~ log1p(nsurv), nsp, axes = F, xlab = "survival", ylab = "death", col = as.numeric(!nsp$trymodel)+1, pch = 17, cex = 0.8)
box()
axis(1, at = log1p(c(0, 1, 10, 100, 1000)), labels = c(0, 1, 10, 100, 1000))
axis(2, at = log1p(c(0, 1, 10, 100, 1000)), labels = c(0, 1, 10, 100, 1000))
abline(h = log1p(x), v = log1p(x))

table(nsp$trymodel)





# Function for fitting models ---------------------------------------------


model_fit = function(data, speciesinfo) {
  
  # create new factor with correct factor levels per species (otherwise problem with margins)
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, " + s(census, bs = 're')", "") 
  form =  as.formula(paste0("mort_next ~ s(dbh, k = k1) + s(all_BA, k = k2)  + s(con_BA, k = k3)"
                            , term_c))
  
  # chose penalty
  # set to default 10 (the same as -1)
  k1 = k2 = 10 
  if (k1 > speciesinfo$unique_dbh) k1 = speciesinfo$unique_dbh - 2
  if (k2 > speciesinfo$unique_all_BA) k2 = speciesinfo$unique_all_BA - 2
  
  # less flexible k for conspecific density
  k3 = 10
  if (k3 > speciesinfo$unique_con_BA) k3 = speciesinfo$unique_con_BA - 2
  
  
  # fit model
  # https://stats.stackexchange.com/questions/27586/spline-df-selection-in-a-general-additive-poisson-model-problem/71300#71300
  mod = try(gam(form
                , family = binomial(link=cloglog)
                , offset = log(interval)
                , data = data
                # , method = "GCV.Cp" # tends to undersmooth
                , method = "REML"
  ) , silent = T
  )
  
  return(mod)
  
  # # Explore model
  # summary(mod)
  # plot(mod, pages = 1, scale = 0, scheme = 2)
  # k.check(mod)
  # 
  # # with gratia
  # draw(mod)
  # draw(derivatives(mod))
  # 
  # # with gamViz
  # vizmod <- getViz(mod)
  # pl = plot(vizmod, nsim = 20) + l_ciLine() + l_fitLine() + l_simLine()
  # print(pl, pages = 1)
  # pl = plot(vizmod, nsim = 20) + l_fitDens() + l_simLine(colour = 1) + theme(legend.position="none")
  # print(pl, pages = 1)
  # 
  # 
  # # Residuals
  # qq(vizmod, rep = 10, method = "auto", CI = "normal", showReps = TRUE,
  #    a.replin = list(alpha = 0.1), discrete = TRUE)
  # 
  # vizmod <- getViz(mod, nsim = 30)
  # gridPrint(check1D(vizmod, "dbh") + l_gridCheck1D(gridFun = sd, showReps = TRUE))
  # gridPrint(check1D(vizmod, "all_BA") + l_gridCheck1D(gridFun = sd, showReps = TRUE))
  # gridPrint(check1D(vizmod, "con_BA") + l_gridCheck1D(gridFun = sd, showReps = TRUE))
  # # red intervals are a 95% posterior credible intervals for the residuals standard deviation, computed using the posterior simulations
  # # the black points are the standard deviation of the observed binned residuals

}


model_convergence = function(model) {
  
  # gam not available
  if (!any(class(model)=="gam")) {
    print(paste(sp, "gam failed"))
  } else {
    
    # gam not converged
    if (!model$converged) {
      print(paste(sp, "no convergence"))
    } else {
      
      # check for complete separation
      # https://stats.stackexchange.com/questions/336424/issue-with-complete-separation-in-logistic-regression-in-r
      # Explore warning "glm.fit: fitted probabilities numerically 0 or 1 occurred"
      eps <- 10 * .Machine$double.eps
      glm0.resids <- augment(x = model) %>%
        mutate(p = 1 / (1 + exp(-.fitted)),
               warning = p > 1-eps,
               influence = order(.hat, decreasing = T))
      infl_limit = round(nrow(glm0.resids)/10, 0)
      # check if none of the warnings is among the 10% most influential observations, than it is okay..
      num = any(glm0.resids$warning & glm0.resids$influence < infl_limit)
      
      # complete separation
      if (num) {
        print(paste(sp, "complete separation is likely"))
      } else {
        
        # missing Vc
        if (is.null(model$Vc)) {
          print(paste(sp, "Vc not available"))
        } else {
        
        # successful model
        return(model)
        }
      }
    }
  }
}







# Fit models --------------------------------------------------------------


res_mod = list()


# Fit models for individual species
for (sp in nsp$sp[nsp$trymodel]) {
  
  # select data for individual species
  dat_sp = dat_mort[dat_mort$sp == sp, ]
  
  # model fit
  mod = model_fit(data = dat_sp, speciesinfo = nsp[nsp$sp == sp, ])
  
  # check model success
  res = model_convergence(model = mod)
  
  # save result
  if (is.character(res)) {
    nsp$rare[nsp$sp == sp] = T  
  } else {
    res_mod[[sp]] = res
  }
}



# Fit models for rare species groups --------------------------------------


# overview rare versus rare and no convergence
table(rare = nsp$rare, trymodel = nsp$trymodel)

# make Rare dependent on stature
nsp$stature = stature$stature[match(nsp$sp, stature$sp)]
nsp$rare_stature = NA
nsp$rare_stature[nsp$rare] = paste0("Rare_", nsp$stature[nsp$rare])


# explore available data per rare species group
dat_mort %>% 
  mutate(sp = nsp$rare_stature[match(sp, nsp$sp)]) %>% 
  filter(!is.na(sp)) %>% 
  group_by(sp) %>% 
  summarise(ndead = sum(surv_next==0),
            nsurv = sum(surv_next==1),
            range_con_BA = max(con_BA) - min(con_BA),
            max_con_BA = max(con_BA),
            unique_con_BA = length(unique(con_BA)),
            unique_all_BA = length(unique(all_BA)),
            unique_dbh = length(unique(dbh))
  ) %>% 
  
  # add issues
  mutate(issue_nval = unique_con_BA < nval,                             # less than nval unique values in conspecific densities
         issue_range = range_con_BA < minrange,                         # range should at least be equal to minrange
         issue_cutoff = ndead < x | nsurv < x,                          # ndead or nsurv below cutoff
         trymodel = !(issue_nval | issue_range | issue_cutoff),         # should a group specific model be tried at all?
  ) -> nsp_rare



# Fit models for species groups
for (sp in unique(nsp_rare$sp[nsp_rare$trymodel])) {
  
  # select data for all species in group
  sps = nsp$sp[which(nsp$rare_stature == sp)]
  dat_sp = dat_mort[dat_mort$sp %in% sps, ]
  
  # model fit
  mod = model_fit(data = dat_sp, speciesinfo = nsp_rare[nsp_rare$sp == sp, ])
  
  
  # check model success
  res = model_convergence(model = mod)
  
  # save result
  if (!is.character(res)) {
    res_mod[[sp]] = res
  }
}





# Regression table via broom::tidy() --------------------------------------

coefs = lapply(res_mod, broom::tidy)
coefs = Map(cbind, coefs, sp = names(coefs))
coefs = do.call(rbind, coefs)



# Model summary via broom::glance() ---------------------------------------

sums = lapply(res_mod, broom::glance)
sums = Map(cbind, sums, sp = names(sums))
sums = do.call(rbind, sums)

# AUC
aucs = lapply(res_mod, function(x) {
  roc <- performance::performance_roc(x, new_data = x$model)
  bayestestR::area_under_curve(roc$Spec, roc$Sens)
})
sums$AUC = unlist(aucs)



# Settings for AMEs -------------------------------------------------------


# fixed value for interval (offset)
# different change settings for conspecific densities
# averaged over all other predictors (total density, dbh...)

interval = 1
change = list(equilibrium = data.frame(con_BA = "paste('+', additive)")
              , invasion = data.frame(con_BA = "c(0, additive)"))
iter = 500





# Absolute AMEs -----------------------------------------------------------


# Calculate absolute AMEs based on manual function get_AME
AME = data.frame()
AMEsamples = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][, i]))
              , iterations = iter
              , samples = T
      )
    }
    )
    
    # AME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change = j, sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    AME = rbind(AME, tempAME)
    
    # AME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    AMEsamples = rbind(AMEsamples, tempSamples)
  }
}



# Relative AMEs -----------------------------------------------------------


# Calculate relative AMEs based on manual function get_AME
rAME = data.frame()
rAMEsamples = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][, i]))
              , iterations = iter
              , relative = T
              , samples = T
      )
    }
    )
    
    # rAME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change = j, sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    rAME = rbind(rAME, tempAME)
    
    # rAME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    rAMEsamples = rbind(rAMEsamples, tempSamples)
  }
}




# plot splines in pdf -----------------------------------------------------

pdf(paste0(path_mortality, run, "/", site, "_mortality.pdf"))
for (i in 1:length(res_mod))  {
  
  # with mgcv
  # plot(res_mod[[i]], pages = 1, scale = 0, main = names(res_mod)[i], all.terms = T) 
  
  # with gamViz
  vizmod <- getViz(res_mod[[i]], post = T, unconditional = T)
  pl = plot(vizmod, nsim = 20, allTerms = T) + 
    l_ciLine() + l_fitLine() + l_simLine() + 
    l_ciBar() + l_fitPoints(size = 1) + 
    l_rug() + 
    # l_points() +
    labs(title = names(res_mod)[i]) 
  print(pl, pages = 1)
  
}
dev.off()







# plot residuals ----------------------------------------------------------



sums$test.spatial = NA
sums$MoransIobs = NA
sums$MoransIexp = NA
sums$MoransIsd = NA


pdf(paste0(path_mortality, run, "/", site, "_residuals.pdf"), width = 10)

for (i in 1:length(res_mod))  {

  sp = names(res_mod)[i]

  # get data (for coordinates)
  dat_sp = dat_mort[dat_mort$sp == sp, ]
  if (nrow(dat_sp) == 0) {   # for rare species groups
    sps = nsp$sp[which(nsp$rare_stature == sp)]
    dat_sp = dat_mort[dat_mort$sp %in% sps, ]
  }

  # get model
  mod = res_mod[[sp]]
  
  # get simulated residuals with DHARMa
  res = simulateResiduals(mod)
  
  # Classic residual checks
  par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
  plotResiduals(res, quantreg = T)
  plotResiduals(res, form = mod$model$dbh, quantreg = T, rank = F, xlab = "dbh (mm)")
  plotResiduals(res, form = mod$model$all_BA, quantreg = T, rank = T, xlab = "all_BA (rank transformed)")
  plotResiduals(res, form = mod$model$con_BA, quantreg = T, rank = T, xlab = "con_BA (rank transformed)")
  plotResiduals(res, form = exp(mod$model$'(offset)'), quantreg = T, rank = F, xlab = "interval (yrs)")

  # recalculate residuals for unique coordinates
  res_spatial = recalculateResiduals(res, paste(dat_sp$gx, dat_sp$gy))
  group_coor = do.call("rbind", strsplit(unique(res_spatial$group), split = " "))

  # Spatial residual checks
  testSpatial = testSpatialAutocorrelation(simulationOutput = res_spatial
                                           , x = group_coor[, 1]
                                           , y = group_coor[, 2]
                                           , plot = F)

  # assign results to sums
  sums$test.spatial[sums$sp == sp] = testSpatial$p.value
  sums$MoransIobs[sums$sp == sp] = testSpatial$statistic["observed"]
  sums$MoransIexp[sums$sp == sp] = testSpatial$statistic["expected"]
  sums$MoransIsd[sums$sp == sp] = testSpatial$statistic["sd"]

  col = colorRamp(c("red", "white", "blue"))(res_spatial$scaledResiduals)
  plot(group_coor[, 1], group_coor[, 2], col = rgb(col, maxColorValue = 255)
       , main = testSpatial$method, cex.main = 0.8
       , xlab = "x", ylab = "y", asp = 1, cex = 0.5
  )

  mtext(sp, outer = T)

}
dev.off()






# Save results ------------------------------------------------------------


save(list = c("AME", "AMEsamples", "rAME", "rAMEsamples", "nsp", "nsp_rare", "coefs", "sums")
     , file = paste0(path_mortality, run, "/", site, "_mortality.Rdata"))

rm(list = ls()[!(grepl("site", ls()) | grepl("path", ls()) | ls()=="run" | ls()=="t0" )])



