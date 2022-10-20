

# Mortality analyses
# Species-specific models with decay optimization: BA





# Packages and data -------------------------------------------------------

library(effects)
library(dplyr)
library(tidyr)
library(broom)
library(mgcv)
library(mgcViz)
library(performance)
library(bayestestR)


source("code/functions_marginal_effects_gam.R")


# Settings ----------------------------------------------------------------


# chose decay type
decay_type = "expn"

dec_fun <- function(sigma, distance, type) {  # define decay function explicitly
  
  if (type == "exp") {
    return( exp(-(1/sigma)*distance) )  
  }
  
  if (type == "expn") {
    return( exp(-(1/(sigma^2))*(distance^2)) )  
  }
}


# chose predictors for local density
predictors = c(con_BA = paste0("ba_con_", decay_type, "_", decay_con),
               all_BA = paste0("ba_all_", decay_type, "_", decay_tot))





# Load data ---------------------------------------------------------------


load(paste0("data_prep/data_3a_mortality/", site, "_tree_3a_mortality_", decay_type, ".Rdata"))


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
load(paste0("../ForestGEO_datacleaning@git/data_species/", site, "_species.Rdata"))
dat_mort = dat_mort[dat_mort$sp %in% species$sp[species$fern.palm == "FALSE"], ]

# remove NA and 0 in interval
dat_mort = dat_mort[!is.na(dat_mort$interval) & !(dat_mort$interval == 0), ]

# limit to small trees
dat_mort = dat_mort[dat_mort$dbh < 100, ]

# save census as factor
dat_mort$census = as.factor(dat_mort$census)
table(dat_mort$census)




# Species handling --------------------------------------------------------


# aggregate species information to set k
dat_mort %>% 
  group_by(sp) %>% 
  summarise(ndead = sum(surv_next==0),
            nsurv = sum(surv_next==1),
            range_con_BA = max(con_BA) - min(con_BA),
            max_con_BA = max(con_BA),
            unique_con_BA = length(unique(con_BA)),
            unique_all_BA = length(unique(all_BA)),
            unique_dbh = length(unique(dbh))
  ) -> nsp





# Function for fitting models ---------------------------------------------


model_fit = function(data, speciesinfo) {
  
  # create new factor with correct factor levels per species (otherwise problem with margins)
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, " + s(census, bs = 're')", "")  # elevation term
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
          
          # summary doensn't work (error message of mgcv is "at least one element of lb must be non-zero")
          # occured with expn edo 3.5 22
          test = tryCatch({summary(model); F}, error=function(e) T)
          if (test) {
            print(paste(sp, "summary(model) not available"))
          } else {
            
            # successful model
            return(model)
          }
        }
      }
    }
  }
}






# Fit models --------------------------------------------------------------


res_mod = list()


# Fit models for individual species
for (sp in nsp$sp) {
  
  # select data for individual species
  dat_sp = dat_mort[dat_mort$sp == sp, ]
  
  # model fit
  mod = model_fit(data = dat_sp, speciesinfo = nsp[nsp$sp == sp, ])
  
  # check model success
  res = model_convergence(model = mod)
  
  # save result
  if (is.character(res)) {
    nsp$convergence[nsp$sp == sp] = F
  } else {
    nsp$convergence[nsp$sp == sp] = T 
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

interval = 5
additive = round(pi*((20/1000)/2)^2, 5) # basal area (m2) of one more neighbor with dbh = 20 mm
change = list(equilibrium = data.frame(con_BA = "paste('+', additive)"))







# Absolute AMEs -----------------------------------------------------------


# Calculate absolute AMEs based on manual function get_AME
AME = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][, i]))
      )
    }
    )
    
    # AME
    tempAME = lapply(temp, function(x) x)
    tempAME = Map(cbind, tempAME, change = j, sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    AME = rbind(AME, tempAME)
    
  }
}



# Relative AMEs -----------------------------------------------------------


# Calculate relative AMEs based on manual function get_AME
rAME = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][, i]))
              , relative = T
      )
    }
    )
    
    # rAME
    tempAME = lapply(temp, function(x) x)
    tempAME = Map(cbind, tempAME, change = j, sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    rAME = rbind(rAME, tempAME)
    
  }
}








# Save results ------------------------------------------------------------

save(list = c("AME", "rAME", "nsp", "coefs", "sums")
     , file = paste0("out/mortality_models/"
                     , run, "/"
                     , format(decay_con, nsmall = 1), "_"
                     , sprintf("%02d", decay_tot), "/"
                     , site, "_mortality.Rdata"))

rm(list = ls()[!(grepl("site", ls()) | ls()=="run" | ls()=="t0")])



