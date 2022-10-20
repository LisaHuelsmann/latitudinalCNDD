

# Functions to calculate average marginal effects for mgcv::gam()


library(MASS)


# Useful functions --------------------------------------------------------


# set step on x for numerical derivative
setstep = function(x) {
  eps = .Machine$double.eps
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}





# Average marginal effects ------------------------------------------------

get_AME = function(mod, data, term
                   , change = NULL
                   , at = NULL
                   , offset = 1
                   , relative = F
                   , iterations = 1000
                   , seed = 10
                   , samples = F) {
  
  # Two data frames with shift in term of interest:
  d0 = d1 = data
  
  # if change is NULL, use numerical derivative
  if (is.null(change)) {
    
    d0[[term]] = d0[[term]] - setstep(d0[[term]])
    d1[[term]] = d1[[term]] + setstep(d1[[term]])
    
  } 
  
  # if change includes "+", use explicit additive change
  if (grepl("\\+", paste(change, collapse = "_"))) {
    
    d1[[term]] = d1[[term]] + as.numeric(gsub("\\+", "", change))
    
  } 
  
  # if change has two values, use explicit change
  if (length(change) == 2) {
    
    d0[[term]] = as.numeric(change[1])
    d1[[term]] = as.numeric(change[2])
    
  }
  
  # Fix at values:
  if (!is.null(at)) {
    for (i in names(at))
      d0[[i]] = at[[i]]
      d1[[i]] = at[[i]]
  }
  
  # Matrices for prediction, map coefs to fitted curves
  Xp0 <- predict(mod, newdata = d0, type="lpmatrix")
  Xp1 <- predict(mod, newdata = d1, type="lpmatrix")
  
  # Model settings
  ilink <- family(mod)$linkinv
  beta <- coef(mod)
  # vc <- mod$Vp # Bayesian
  vc <- mod$Vc # Bayesian accounting for smoothing parameter uncertainty
  

  # marginal effects
  pred0   <- 1 - (1-ilink(Xp0 %*% beta))^offset
  pred1   <- 1 - (1-ilink(Xp1 %*% beta))^offset
  ME <- (pred1-pred0)
  
  # if change is NULL, use numerical derivative
  if (is.null(change)) {
    ME <- ME/(d1[[term]] - d0[[term]])
  } 

  
  # convert to relative if requested
  if (relative == T) ME = ME/pred0
  
  # # Alternative for relative:
  # ME1 = ME/pred0; ME2 = ME/mean(pred0)
  # plot(ME, ME2)
  # points(ME, ME1, col = 2)
  # par(mfrow = c(1, 2))
  # hist(ME1); hist(ME2)
  # par(mfrow = c(1, 1))
  
  # average marginal effect
  AME = mean(ME)
  
  
  # generate AME samples
  
  # variance of average marginal effect via "posterior" simulation
  # simulate from multivariate normal using model beta means and covariance matrix
  if (!is.null(seed)) set.seed(seed)
  coefmat = mvrnorm(n = iterations
                    , mu = beta
                    , Sigma = vc)
  
  # estimate AME from each simulated coefficient vector
  AMEs = apply(coefmat, 1, function(coefrow) {
    
    # marginal effects
    pred0   <- 1 - (1-ilink(Xp0 %*% coefrow))^offset
    pred1   <- 1 - (1-ilink(Xp1 %*% coefrow))^offset
    ME <- (pred1-pred0)
    
    # if change is NULL, use numerical derivative
    if (is.null(change)) {
      ME <- ME/(d1[[term]] - d0[[term]])
    } 
    
    # convert to relative if requested
    if (relative == T) ME = ME/pred0
    
    # average marginal effect
    AME = mean(ME)
    return(AME)
  })
  
  # Combine results
  if (!samples) {
    res = data.frame(term
                     , estimate = AME
                     , std.error = sqrt(var(AMEs))  # not a good idea
                     , estimate.sim = mean(AMEs)    # not a good idea
                     , offset
                     , change.value = paste(change, collapse = "_"))
    return(res) 
    
  } else {
    
    res_sums = data.frame(term
                     , estimate = AME
                     , std.error = sqrt(var(AMEs)) # not a good idea
                     , offset
                     , change.value = paste(change, collapse = "_"))
    
    res_samples = data.frame(term
                             , estimate = AMEs
                             , MLE = AME
                             , offset
                             , change.value = paste(change, collapse = "_"))
    res = list(res_sums, res_samples)
    return(res)  
    
  }
}










