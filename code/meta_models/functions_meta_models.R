


# Functions to analyze outputs






# Aesthetic functions -----------------------------------------------------


# function to an alpha value to a colour
add.alpha <- function(col=NULL, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# https://github.com/cran/geometry/blob/master/R/polyarea.R
polyarea <- function(x, y, d=1) {
  if (is.vector(x) & is.vector(y)) {
    if (length(x) == length(y)) {
      a <- abs(sum(x*(magic::shift(y, -1) - magic::shift(y, 1))))/2
    } else {
      stop("x and y must have the same shape")
    }
  } else {
    if (is.array(x) & is.array(y)) {
      if (length(dim(x)) != 2) {
        stop("Arrays must have two dimensions.")
      }
      if (all(dim(x) == dim(y))) {
        v <- c(0, 0)
        v[d] <- 1
        a <- abs(apply(x*(magic::ashift(y, -v) - magic::ashift(y, v)), 3 - d, sum))/2
      } else {
        stop("x and y must have the same shape")
      }
    } else {
      stop("x and y must be of same type")
    }
  }
  names(a) <- NULL
  return(a)
}


# coordinates for plots
circular_coor = function(order, inner) {
  
  # select corners
  corners = 1:4*round(length(order)/4, 0)-1
  cornerorder = order[corners]
  cornerorder = cornerorder[!is.na(cornerorder)]
  
  # reduce order for circular plotting
  circorder = order[-corners]
  
  # settings
  deg2rad <- function(deg) {(deg * pi) / (180)}
  twist = 63
  sizeMiddle = 0.7
  sideMid = 0.8*(1-sizeMiddle)
  sideDis = 0.3*(1-sizeMiddle)
  
  # object for coordinates
  circCoords = matrix(0, length(circorder) + 1, 4)
  circCoords[1,] =  mid = rep(c(0.5 - 0.5*sizeMiddle, 0.5 + 0.5*sizeMiddle), 2)
  
  # inner circle
  angles_inner = seq(0+twist+10,360+twist+10,length.out = inner+1)[1:inner]
  for(p in 2:(inner+1)){
    pp = p - 1
    circCoords[p,1] = cos(deg2rad(angles_inner[pp])) * (sideMid - 0.01) + 0.5 - sideDis
    circCoords[p,2] = circCoords[p,1] + 2*sideDis
    if(circCoords[p,1] > circCoords[p,2]) circCoords[p,] = circCoords[p,c(2,1,3,4)]
    circCoords[p,3] = sin(deg2rad(angles_inner[pp])) * (sideMid - 0.01) + 0.5 -sideDis
    circCoords[p,4] = circCoords[p,3] + 2*sideDis
    if( circCoords[p,4] > 1) circCoords[p,4] = 1
  }
  
  # outer circle
  angles_outer = seq(0+twist,360+twist,length.out = length(circorder)-inner+1)[1:(length(circorder) - inner)]
  for(p in (inner+2):(length(circorder)+1)){
    pp = p - (inner+1)
    circCoords[p,1] = cos(deg2rad(angles_outer[pp])) * (sideMid + 0.17) + 0.5 - sideDis
    circCoords[p,2] = circCoords[p,1] + 2*sideDis
    if(circCoords[p,1] > circCoords[p,2]) circCoords[p,] = circCoords[p,c(2,1,3,4)]
    circCoords[p,3] = sin(deg2rad(angles_outer[pp])) * (sideMid + 0.17) + 0.5 -sideDis
    circCoords[p,4] = circCoords[p,3] + 2*sideDis
    if( circCoords[p,4] > 1) circCoords[p,4] = 1
  }
  
  # bring in circular order
  circCoords = circCoords[order(c(angles_inner, angles_outer))+1, ]
  
  # add corners
  xcorners = c(0.1, 0.1, 0.9, 0.9)[1:length(cornerorder)]
  ycorners = c(0.9, 0.1, 0.1, 0.9)[1:length(cornerorder)]
  cornerCoords = matrix(0, length(cornerorder), 4)
  cornerCoords[, 1] = xcorners-sideDis
  cornerCoords[, 2] = xcorners+sideDis
  cornerCoords[, 3] = ycorners-sideDis
  cornerCoords[, 4] = ycorners+sideDis
  
  # combine coordinates
  multiCoords = matrix(0, length(order), 4)
  multiCoords[-corners, ] = circCoords
  multiCoords[corners[1:length(cornerorder)], ] = cornerCoords
  
  return(multiCoords)
  
}


format_pval = function(p, addp = T) {
  
  p_return = sfsmisc::pretty10exp(p, sub10 = c(-4, 4), digits = 2)
  
  if (addp) {
    for (i in seq(along = p_return)) {
      p_return[[i]] <- substitute(p==pvalue, list(pvalue = p_return[[i]]))
    }
  }
  
  return(p_return)
}


decimalplaces <- function(x) {
    nchar(strsplit(sub('0+$', '', x), ".", fixed = TRUE)[[1]][[2]])
}




# Fit meta-regressions ----------------------------------------------------


# Site specific models without and with moderator
# includes tryCatch to record errors and likelihood profile problems

rma_fit = function(dat, moderator = ~ 1) {
  
  # utility function
  f = function(dat, moderator) {
    rma(yi = yi
        , vi = vi
        , mods = moderator
        , method = "REML"
        , data = dat
    )
  }
  
  mod = tryCatch(f(dat = dat, moderator = moderator)
                 , warning = function(w) return(list(res = f(dat = dat, moderator = moderator)
                                                     , warning = paste(w, "okay")))
                 , error = function(e) e)
  
  return(mod)
  
}


# model convergence and other issues
rma_convergence = function(mod) {
  
  # no warnings or errors
  if(any(class(mod) == "rma")) {
    reliable = T
  } else {
    
    # error
    if (any(class(mod)=="error")) {
      reliable = F
    } else {
      
      # warning for likelihood profile
      if (grepl("profile", mod$warning)) {
        reliable = F
      } else {                 
        reliable = T
        mod = mod$res
      }
    }   
  }
  
  return(list(mod = mod, reliable = reliable))
  
}







# Plot CNDD estimates  ----------------------------------------------------


# plot estimated effects against abundance
plot_estimates_vsAbundance = function(data, type, order, names = NULL
                                      , trans, backtrans
                                      , multiCoords
                                      , cols = F, color.axes = "black"
                                      , cex.text = 0.4
                                      , markRare = F
                                      , returnMod = F
                                      ) {
  
  # limits
  xlim = log(c(min(data$abundance), max(data$abundance)))
  ylim = quantile(data$estimate, probs = c(0.2, 0.8))
  
  # y sequence
  digitrounding = ifelse(type == "rAME", 100, 1000)
  yspacing = 2
  yseq = floor(digitrounding*backtrans(ylim))/digitrounding       # backtransformed 
  yseq = c(yseq[rep(1, yspacing)], 0,  yseq[rep(2, yspacing)])
  yseq = yseq/c(1, yspacing, 1, yspacing, 1)
  
  
  
  # use site names from order when not provided in names
  if (is.null(names)) names = order
  
  
  if (returnMod) models = list()
  
  for (site in order) {
    
    i = which(order == site)
    
    
    if (!any(data$site == site)) {
      
      plot.new()
      
    } else {
      
      
      
      # prepare escalc for site
      dat_meta = escalc(measure = "GEN"
                        , yi = estimate
                        , sei = std.error
                        , slab = sp
                        , data = data[data$site==site, ]
      )
      
      
      # fit model
      mod = rma_fit(dat_meta, moderator = ~ log(abundance))
      
      # check for non-reliable results 
      check = rma_convergence(mod)
      mod = check$mod
      reliable = check$reliable
      
      
      
      # generate predictions and save model when model okay
      if (reliable) {
        npred = 500
        xs <- seq(min(log(dat_meta$abundance)), max(log(dat_meta$abundance)), length = npred)
        sav <- predict(mod, newmods = t(xs))
        
        # save model if returnMod
        if (returnMod) models[[site]] = mod
      }
      
      # significance
      significant = dat_meta$significant
      
      # plot CNDD against abundance
      par(fig = multiCoords[which(order == site),], new = T)
      plot(-100, -100
           , xlim = xlim
           , ylim = ylim
           , xlab = ""
           , ylab = ""
           , xaxt = "n"
           , yaxt = "n"
           , bty='n'
      )
      
      # white background
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = "white", border = NA)
      
      # axes
      abline(h = 0, col = "lightgrey", lwd = 0.6, lty = 2)
      axis(1, log(c(0.01, 0.1, 1, 10, 100, 1000, 10000)), labels = NA
           , tck = -0.025, col = color.axes, lwd = 0.6, lwd.ticks = 0.6)
      axis(2, trans(yseq), labels = NA
           , tck = -0.025, col = color.axes, lwd = 0.6, lwd.ticks = 0.6)
      
      # axes labels
      axis(1, log(c(0.1, 10, 1000)), labels = c(0.1, 10, 1000), tick = F,
           cex.axis = cex.text, lwd = 0, line = -1.4, col.axis = color.axes)
      axis(2, trans(yseq), labels = 100*yseq, las = 1, lwd = 0
           , cex.axis = cex.text, line = -0.8, col.axis = color.axes)  
      
      
      # add points of individual estimates
      if (markRare) pch = ifelse(grepl("Rare", dat_meta$sp), 18, 16)
      if (!markRare) pch = 16
      points(log(dat_meta$abundance), dat_meta$yi
             , pch = pch
             , cex = 0.18*(scale(1/(dat_meta$vi)^(1/10)) + min(1/(dat_meta$vi)^(1/10)))
             , col= rgb(0, 0, 0, 0.4*significant + 0.2))
      
      # add polygon and line for model fit
      if (reliable) {
        confheight = mean(sav$ci.ub-sav$ci.lb)
        plotheight = diff(ylim)
        
        alpha = 1-8*(confheight/plotheight)
        if (alpha < 0.2) alpha = 0.2
        
        polygon(c(xs, rev(xs)), c(sav$ci.lb, rev(sav$ci.ub))
                , col = add.alpha(cols[i], alpha)
                , border = NA)
        lines(xs, sav$pred, lwd=1)
        
        # p-value
        text(min(xs), sav$ci.ub[1]+0.02*ylim[2]
             , format_pval(summary(mod)$pval[2])
             , cex = cex.text, adj = c(0, 0)
        )
      }
      
      # add how much of the data (in terms of nobs) is outside plotting area
      outside = dat_meta$estimate < ylim[1] | dat_meta$estimate > ylim[2]
      if (any(outside)) {
        dat_outside = sums_global$nobs[sums_global$site == site & sums_global$sp %in% dat_meta$sp[outside]]
        dat_outside = sum(dat_outside)/sum(sums_global$nobs[sums_global$site == site])
        text(xlim[1], ylim[1]+0.02*ylim[2], paste0(round(100*dat_outside, 0), "% of data outside")
             , cex = cex.text, adj = c(0, 0), col = "grey70")
      } else {
        text(xlim[1], ylim[1]+0.02*ylim[2], "no data outside"
             , cex = cex.text, adj = c(0, 0), col = "grey70")
      }
      par(xpd = NA, new = T)
      text(0.95*xlim[1], ylim[2]-0.3*ylim[2]
           , gsub('(.{1,30})(\\s|$)', '\\1\n', names[i])
           , cex = cex.text, adj = c(0, 0), font = 2, lheight=.6)
      par(xpd = F, new = T)
    }
  }
  if (returnMod) return(models)
}




# Get predictions -----------------------------------------------------


# return predictions for object x
# either against latitude or against abundance (two separate function)
# x is a list that contains models, data and transformations for one dataset and CNDD definition


get_predictions_latitude = function(x
                                    , select = "mod0"
                                    , abundances = c(1, 10, 1000)
                                    , pvalue = F) {
  
  # extract model
  model = x[[select]]
  
  backtrans = get(paste0("backtrans_", x$type))
  
  # global y limits
  y = 0.8*100*c(predict(model, transf = backtrans)$ci.lb, 
            predict(model, transf = backtrans)$ci.ub)
  ylim = range(y)
  
  res = list()
  
  # produce predictions and other necessary info for plotting per abundance
  for (i in 1:length(abundances)) {
    
    # xlimits, xvalues and new data
    xlim = c(0, max(x$data$latitude))
    x_line = seq(xlim[1], xlim[2], length.out = 500)
    x_poly = c(x_line, rev(x_line))
    newDat = cbind(x$transLat(x_line, x$ref_lat)
                   , x$transAbund(abundances[i], x$ref_abund)
                   , x$transLat(x_line, x$ref_lat)*x$transAbund(abundances[i], x$ref_abund)
    )
    newDat = newDat[, 1:(length(model$coef)-1)] # remove interaction column when not needed
    
    # predictions 
    predictions = predict(model, newmods = newDat, transf = backtrans)
    pred_line = 100 * predictions$pred
    pred_poly = 100 * c(predictions$ci.lb, rev(predictions$ci.ub))
    
    # p-values (may require model refitting)
    if (pvalue) {
      if (abundances[i] == x$ref_abund) {
        p_value = model$pval[rownames(model$b)=="tLatitude"]
      } else {
        
        # refit model with new ref_abund
        refit_data = x$data
        refit_data$tAbundance =  x$transAbund(refit_data$abundance, ref_abund = abundances[i])
        refit_model = metafor::update.rma(model, data = refit_data)
        p_value = refit_model$pval[rownames(refit_model$b)=="tLatitude"]
      }
      p_value = format_pval(p_value)
      
    } else {
      p_value = NULL
    }
    
    # compile output
    res[[i]] = list(xlim = xlim, ylim = ylim
                    , pred_line = pred_line, x_line = x_line
                    , pred_poly = pred_poly, x_poly = x_poly
                    , pvalue = p_value
                    , type = x$type
                    , label = paste("CNDD assessed as", x$type, "in", x$term, "calculated at", 
                                    x$change, "densities."))
  }
  names(res) = abundances
  return(res)
}

get_predictions_abundance = function (x
                                      , select = "mod0"
                                      , latitudes = c(0, 25, 50)
                                      , pvalue = F) {
  
  # extract model
  model = x[[select]]
  
  backtrans = get(paste0("backtrans_", x$type))
  
  # global y limits
  y = 0.8*100*c(predict(model, transf = backtrans)$ci.lb, 
            predict(model, transf = backtrans)$ci.ub)
  ylim = range(y)
  
  res = list()
  
  # produce predictions and other necessary info for plotting per latitude
  for (i in 1:length(latitudes)) {
    
    # xlimits, xvalues and new data
    xlim = c(min(x$data$abundance), max(x$data$abundance))
    x_line = seq(xlim[1], xlim[2], length.out = 250)
    x_poly = c(x_line, rev(x_line))
    newDat = cbind(x$transLat(latitudes[i], x$ref_lat)
                   , x$transAbund(x_line, x$ref_abund)
                   , x$transLat(latitudes[i], x$ref_lat) * x$transAbund(x_line, x$ref_abund))
    
    # predictions 
    predictions = predict(model, newmods = newDat, transf = backtrans)
    pred_line = 100 * predictions$pred
    pred_poly = 100 * c(predictions$ci.lb, rev(predictions$ci.ub))
    
    # p-values (may require model refitting)
    if (pvalue) {
      if (latitudes[i] == x$ref_lat) {
        p_value = model$pval[rownames(model$b)=="tAbundance"]
      } else {
        
        # refit model with new ref_abund
        refit_data = x$data
        refit_data$tLatitude =  x$transLat(refit_data$latitude, ref_lat = latitudes[i])
        refit_model = metafor::update.rma(model, data = refit_data)
        p_value = refit_model$pval[rownames(refit_model$b)=="tAbundance"]
      }
      p_value = format_pval(p_value)
      
    } else {
      p_value = NULL
    }
    
    # compile output
    res[[i]] = list(xlim = xlim, ylim = ylim
                    , pred_line = pred_line, x_line = x_line
                    , pred_poly = pred_poly, x_poly = x_poly
                    , pvalue = p_value
                    , type = x$type
                    , label = paste("CNDD assessed as", x$type, "in", x$term, "calculated at", 
                                    x$change, "densities."))
  }
  names(res) = latitudes
  return(res)
}




# Plot CNDD predictions ---------------------------------------------------


# based on previously generated predictions
# either against latitude or against abundance (two separate function)


# against latitude for defined abundances
plot_latitude = function(preds
                         , names = NULL
                         , col_run = NULL
                         , labelsx = "outer"
                         , labelsy = "outer"
                         , xlim = NULL
                         , ylim = NULL
                         , panel = NULL) {
  
  runs = names(preds)
  abundances = as.numeric(names(preds[[1]]))
  type = preds[[1]][[1]]$type
  
  # generate color if needed
  if (is.null(col_run)) col = "orchid4"
  
  # combined x and y range
  if (is.null(xlim)) xlim = range(unlist(lapply(preds, function(x) x[[1]]$xlim)))
  if (is.null(ylim)) ylim = range(unlist(lapply(preds, function(x) x[[1]]$ylim)))
  
  
  # abundances
  for (i in 1:length(abundances)) {
    
    # plotting
    plot(-10, -10, xlim = xlim, ylim = ylim
         , xlab = "", ylab = ""
         , main = ""
         , bty='l'
         , xaxt = "n"
         , yaxt = "n"
    )
    
    # add label for panel if needed
    if (i == 1 & !is.null(panel)) text(grconvertX(0.02, from = "ndc")
                                       , grconvertY(0.98, from = "ndc")
                                       , labels = panel
                                       , cex = 3/2
                                       , xpd=NA)
    
    # abundance label
    if (length(abundances) > 1) {
      graphics::text(xlim[2], 0.9*ylim[2]
                     , paste0("Species abundance\n", abundances[i], " N/ha")
                     , adj = 1, pos = 2)
    }
    
    # axes
    axis(1)
    if (i == 1) axis(2) else axis(2, labels = F)
    abline(h = 0, col = "grey70")
    
    # inner labels
    if (i %in% labelsx) {
      mtext("absolute latitude (°)", 1, 2.5, cex = 2/3)
    }
    if (i %in% labelsy) {
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      mtext(paste("stabilizing CNDD", unit), 2, 3.5, las = 0, cex = 2/3)
    }
    
    for (run in runs) {
      
      # extract prediction
      pred_run = preds[[run]][[as.character(abundances[i])]]
      
      # define color
      if (!is.null(col_run)) col = col_run[which(runs == run)]
      
      lines(pred_run$x_line
            , pred_run$pred_line
            , col = col
      )
      polygon(pred_run$x_poly
              , pred_run$pred_poly
              , col = add.alpha(col, 0.1)
              , border = F)
      
      # p-value
      if (!is.null(pred_run$pvalue)) {
        par(xpd = NA)
        
        # define location
        xloc = pred_run$x_poly[1]
        yloc = pred_run$pred_line[1] + 0.8*(rev(pred_run$pred_poly)[1] - pred_run$pred_line[1])
        yloc = ifelse(yloc < ylim[2], yloc, ylim[2])
        yloc = ifelse(yloc > ylim[1], yloc, ylim[1])
        
        # print p-value
        text(xloc
             , yloc
             , pred_run$pvalue
             , cex = 1
             , adj = c(0, 1))
        par(xpd = F)
      }
      
      # add run name if i == 1
      if (i == 1) {
        par(xpd = NA)
        
        # define location for labels
        xloc = min(pred_run$x_poly)
        yloc = ifelse(which(runs == run) != 3
                      , pred_run$pred_poly[length(pred_run$pred_poly)]
                      , pred_run$pred_poly[1])
        yloc = ifelse(yloc < ylim[2], yloc, ylim[2])
        yloc = ifelse(yloc > ylim[1], yloc, ylim[1])
        text(xloc
             , yloc
             , names[which(runs == run)]
             , col = col
             , adj = c(0, ifelse(which(runs == run) == 3, 1, 0))
             , cex = 0.8
        )
        par(xpd = F)
      }
    }
  }
  
  # labels
  # mtext(preds[[run]][[1]]$label, 3, 1, outer = T)
  
  if (labelsx == "outer") {
    mtext("absolute latitude (°)", 1, 1, outer = T, cex = 2/3)
  }
  if (labelsy == "outer") {
    unit = ifelse(type == "AME", "(% / year)", "(%)")
    mtext(paste("stabilizing CNDD", unit), 2, 2.5, outer = T, las = 0, cex = 2/3)
  }
  
}


# against abundance for defined latitudes
plot_abundance = function(preds
                          , names = NULL
                          , col_run = NULL        # for different runs
                          , col_lat = NULL        # for different latitudes
                          , latitude_names = NULL
                          , labelsx = "outer"
                          , labelsy = "outer"
                          , xlim = NULL
                          , ylim = NULL
                          , panel = NULL) {
  
  runs = names(preds)
  latitudes = as.numeric(names(preds[[1]]))
  type = preds[[1]][[1]]$type
  
  # generate color if needed
  if (is.null(col_run) & is.null(col_lat)) col = "orchid4"
  
  # combined x and y range
  if (is.null(xlim)) xlim = range(unlist(lapply(preds, function(x) x[[1]]$xlim)))
  if (is.null(ylim)) ylim = range(unlist(lapply(preds, function(x) x[[1]]$ylim)))
  
  
  for (i in 1:length(latitudes)) {
    
    # plotting
    plot(-10, -10, xlim = transAbund(xlim, ref_abund), ylim = ylim
         , xlab = "", ylab = ""
         , xaxt = "n"
         , yaxt = ifelse(latitudes[i] == latitudes[1], "s", "n"), xaxt = "n"
         , bty='l'
    )
    
    # add label for panel if needed
    if (i == 1 & !is.null(panel)) text(grconvertX(0.02, from = "ndc")
                                       , grconvertY(0.48, from = "ndc")
                                       , labels = panel
                                       , cex = 3/2
                                       , xpd=NA)
    
    # latitude label
    if (is.null(latitude_names)) {
      graphics::text(transAbund(xlim, ref_abund)[2], 0.9*ylim[2], paste0("Latitude\n", latitudes[i], "°"), 
                     adj = 1, pos = 2)
    } else {
      graphics::text(transAbund(xlim, ref_abund)[2], 0.9*ylim[2], paste0(latitude_names[i], "\n", latitudes[i], "°"), 
                     adj = 1, pos = 2)
    }
    
    # axes
    axis(1, transAbund(c(0.1, 1, 10, 100, 1000, 10000), ref_abund), labels = F)
    axis(1, transAbund(c(0.1, 10, 1000), ref_abund), labels = c(0.1, 10, 1000), tick = F)
    
    if (latitudes[i] == latitudes[1]) axis(2) else axis(2, labels = F)
    abline(h = 0, col = "grey70")
    
    
    # inner labels
    if (i %in% labelsx) {
      mtext("abundance (N/ha)", 1, 2.5, cex = 2/3)
    }
    if (i %in% labelsy) {
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      mtext(paste("stabilizing CNDD", unit), 2, 3.5, las = 0, cex = 2/3)
    }
    
    
    for (run in runs) {
      
      # extract prediction
      pred_run = preds[[run]][[as.character(latitudes[i])]]
      
      # define color
      if (!is.null(col_run)) col = col_run[which(runs == run)]
      if (!is.null(col_lat)) col = col_lat[i]
      
      lines(transAbund(pred_run$x_line, ref_abund)
            , pred_run$pred_line
            , col = col
      )
      polygon(transAbund(pred_run$x_poly, ref_abund)
              , pred_run$pred_poly
              , col = add.alpha(col, 0.1)
              , border = F)
      
      # p-value
      if (!is.null(pred_run$pvalue)) {
        par(xpd = NA)
        
        # define location
        xloc = transAbund(pred_run$x_poly[1], ref_abund)
        yloc = pred_run$pred_line[1] + 0.8*(rev(pred_run$pred_poly)[1] - pred_run$pred_line[1])
        yloc = ifelse(yloc < ylim[2], yloc, ylim[2])
        yloc = ifelse(yloc > ylim[1], yloc, ylim[1])
        
        # print p-value
        text(xloc
             , yloc
             , pred_run$pvalue
             , cex = 1
             , adj = c(0, 1))
        par(xpd = F)
      }
      
      # add run name if i == 1
      if (i == 1) {
        par(xpd = NA)
        
        # define location for labels
        xloc = min(transAbund(pred_run$x_poly, ref_abund))
        yloc = ifelse(which(runs == run) != 3
                      , pred_run$pred_poly[length(pred_run$pred_poly)]
                      , pred_run$pred_poly[1])
        yloc = ifelse(yloc < ylim[2], yloc, ylim[2])
        yloc = ifelse(yloc > ylim[1], yloc, ylim[1])
        text(xloc
             , yloc
             , names[which(runs == run)]
             , col = col
             , adj = c(0, ifelse(which(runs == run) == 3, 1, 0))
             , cex = 0.8
        )
        par(xpd = F)
      }
    }
    
  }
  
  # labels
  mtext(preds[[run]][[1]]$label, 3, 1, outer = T)
  
  if (labelsx == "outer") {
    mtext("abundance (N/ha)", 1, 1, outer = T, cex = 2/3)
  }
  if (labelsy == "outer") {
    unit = ifelse(type == "AME", "(% / year)", "(%)")
    mtext(paste("stabilizing CNDD", unit), 2, 2.5, outer = T, las = 0, cex = 2/3)
  }
  
}


