


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




# Labeling points functions ----------------------------------------------


# originally from https://github.com/JosephCrispell/basicPlotteR/blob/master/R/addNonOverlappingTestLabelsOrPoints.R


# Tutorials
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
#http://r-pkgs.had.co.nz/description.html
#https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html

## Packages to install
#install.packages("devtools")
#install.packages("digest")
#devtools::install_github("klutometis/roxygen")

## Packages to load
#library("devtools")
#library("roxygen2")

## Creating package
#packageDirectory <- "/home/josephcrispell/Desktop/Research/basicPlotteR/"
#usethis::create_package(packageDirectory)

## Documenting changes
#setwd(packageDirectory)
#document()

## Install
#setwd("..")
#install("basicPlotteR")

#' Add non-overlapping text labels to plot
#'
#' This function is similar to the \code{text()} function but it will attempt to re-locate labels that will overlap
#' @param xCoords A vector containing the X coordinates for labels
#' @param yCoords A vector containing the Y coordinates for labels
#' @param labels A vector containing the labels to be plotted
#' @param cex.label A number to scale the size of the plotted labels. Defaults to 1
#' @param col.label The colour of the plotted labels. Defaults to "red". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.line The colour of the line to plot from relocated labels to original location. Defaults to "black". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.background An optional colour for a background polygon plotted behind labels. Defaults to NULL - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param lty A number detailing the type of line to plot from relocated labels to original location. 0: blank, 1: solid, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, and 6: twodash. Defaults to 1. Multiple line types can be provided. If more options than labels provided types will be recycled.
#' @param lwd A number to scale the size of line from relocated labels to original location. Defaults to 1. Multiple line widths can be provided. If more options than labels provided widths will be recycled.
#' @param border The colour of the border to be plotted around the polygon. Defaults to NA - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param avoidPoints A logical variable indicating whether labels shouldn't be plotted on top of points. Defaults to TRUE
#' @param keepLabelsInside A logical variable indicating whether the labels shouldn't be plotted outside of plotting region. Defaults to TRUE
#' @param cex.pt A number used to scale the points plotted on the graph that labels are to be added to. Defaults to 1
#' @keywords text label plot
#' @export
#' @examples
#' # Create some random points
#' n <- 50
#' coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")
#'
#' # Plot them without labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#'
#' # With potentially overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE)
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.label="black")
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.background=rgb(0,0,0, 0.75), col.label="white")
addTextLabels <- function(xCoords, yCoords, labels, cex.label=1, col.label="red", col.line="black", col.background=NULL,
                          lty=1, lwd=1, border=NA, avoidPoints=TRUE, keepLabelsInside=TRUE, cex.pt=1){
  
  
  # Check that the input data are in the correct format #
  
  
  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addTextLabels() The vectors containing the X and Y coodinates must be the same length.")
  }
  if(length(xCoords) != length(labels)){
    stop("addTextLabels() The vector of labels must be the same length as the coordinate vectors.")
  }
  
  
  # Get the axis limits #
  
  
  # Get the axis limits
  axisLimits <- graphics::par("usr")
  
  
  # Check for NA coordinates #
  
  
  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){
    
    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")
    
    # Check for each of the parameters that can have multiple parameters
    if(length(col.line) == length(xCoords)){
      col.line = col.line[-indicesOfNAs]
    }
    if(length(col.background) == length(xCoords)){
      col.background = col.background[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty = lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd = lwd[-indicesOfNAs]
    }
    if(length(border) == length(xCoords)){
      border = border[-indicesOfNAs]
    }
    
    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]
    
    # Remove the respective labels
    labels <- labels[-indicesOfNAs]
  }
  
  
  # Check if axes are logged #
  
  
  # Check X axis
  xAxisLogged <- FALSE
  if(graphics::par("xlog")){
    
    # Note that X axis was logged
    xAxisLogged <- TRUE
    
    # Log the X coordinates
    xCoords <- log10(xCoords)
    
    # Reset the X axis logged flag - fools points and polygon commands below
    graphics::par(xlog=FALSE)
  }
  
  # Check Y axis
  yAxisLogged <- FALSE
  if(graphics::par("ylog")){
    
    # Note that Y axis was logged
    yAxisLogged <- TRUE
    
    # Log the Y coordinates
    yCoords <- log10(yCoords)
    
    # Reset the Y axis logged flag - fools points and polygon commands below
    graphics::par(ylog=FALSE)
  }
  
  #
  # Store the point information #
  #
  
  # Store the input coordinates and labels
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "Labels"=labels, "N"=length(xCoords), "cex"=cex.pt)
  
  # Set the amount to pad onto height and width
  heightPad <- 0.5
  widthPad <- 0.04
  if(is.null(col.background)){
    heightPad <- 0
    widthPad <- 0
  }
  
  # Calculate the label heights and widths
  pointInfo <- calculateLabelHeightsAndWidths(pointInfo=pointInfo, cex=cex.label,
                                              heightPad=heightPad, widthPad=widthPad)
  

  # Produce a list of alternative locations #

  
  # Generate the alternative locations
  alternativeLocations <- generateAlternativeLocations(axisLimits)
  
  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)
  

  # Create a list to store the information about plotted points #

  
  # Initialise the list to store the information about plotted labels
  plottedLabelInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)
  
  
  # Add labels to plot assigning new locations where necessary #
  
  
  # Plot the point label
  for(i in seq_len(pointInfo$N)){
    
    # Set the colours for plotting the label - allows multiple colours and cycling through colours
    labelColour <- setOption(options=col.label, index=i)
    backgroundColour <- setOption(options=col.background, index=i)
    borderColour <- setOption(options=border, index=i)
    
    # Set the line characteristics
    lineColour <- setOption(options=col.line, index=i)
    lineType <- setOption(options=lty, index=i)
    lineWidth <- setOption(options=lwd, index=i)
    
    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    label <- pointInfo$Labels[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]
    
    # Get a new location
    newLocationIndex <- chooseNewLocation(pointInfo, i, alternativeLocations, distances, plottedLabelInfo, axisLimits, keepLabelsInside)
    
    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 && 
       (avoidPoints == TRUE || tooClose(x, y, height, width, plottedLabelInfo) || outsidePlot(x, y, height, width, axisLimits))){
      
      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]
      
      # Add line back to previous location
      addLineBackToOriginalLocation(altX=altX, altY=altY, x=x, y=y, label=label,
                                    cex=cex.label, col=lineColour, lty=lineType, lwd=lineWidth, heightPad=heightPad,
                                    widthPad=widthPad)
      
      # Add label
      addLabel(x=altX, y=altY, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour, heightPad=heightPad, widthPad=widthPad)
      
      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)
      
      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]
      
    }else{
      
      # Add label
      addLabel(x=x, y=y, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour,
               heightPad=heightPad, widthPad=widthPad)
      
      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)
    }
  }
  
  
  # Return axes logged flags to original state - for if person makes any future plots #
  
  
  graphics::par(xlog=xAxisLogged)
  graphics::par(ylog=yAxisLogged)
  
}

#' Add non-overlapping points to plot
#'
#' This function is similar to the \code{points()} function but it will attempt to re-locate points that will overlap
#' @param xCoords A vector containing the X coordinates for labels
#' @param yCoords A vector containing the Y coordinates for labels
#' @param col.line The colour of the line to plot from relocated points to original location. Defaults to "black". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param lty A number detailing the type of line to plot from relocated labels to original location. 0: blank, 1: solid, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, and 6: twodash. Defaults to 1. Multiple line types can be provided. If more options than labels provided types will be recycled.
#' @param lwd A number to scale the size of line from relocated labels to original location. Defaults to 1. Multiple line widths can be provided. If more options than labels provided widths will be recycled.
#' @param keepInside A logical variable indicating whether the points shouldn't be plotted outside of plotting region. Defaults to TRUE
#' @param cex A number used to scale the size of the points plotted. Defaults to 1
#' @param avoidFactor A number that increases (values > 1) or decreases (values < 1) the amount of space alloted to each point. Defaults to 1
#' @param bg A character string that defines the background colour of our point. Defaults to NULL
#' @param ... Arguments to be passed to the \code{points()} function
#' @keywords points x y plot
#' @export
#' @examples
#' # Create some random points
#' n <- 50
#' coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")
#'
#' # Plot points and allow overlapping
#' plot(x=coords$X, y=coords$Y, bty="n", xaxt="n", yaxt="n", cex=3, xlab="X", ylab="Y")
#'
#' # Plot points and avoid overlapping
#' plot(x=NULL, y=NULL, xlim=range(coords$X), ylim=range(coords$Y), bty="n", xaxt="n", yaxt="n", xlab="X", ylab="Y")
#' addPoints(coords$X, coords$Y, cex=3, col.line="red")
addPoints <- function(xCoords, yCoords, col.line="black", lty=1, lwd=1, keepInside=TRUE, cex=1, avoidFactor=1, bg=NULL,
                      ...){
  

  # Create cex/bg value for each point #

  
  # Note the cex is to be applied to each point
  if(length(cex) != length(xCoords)){
    cex <- rep(cex, ceiling(length(xCoords)/length(cex)))[1:length(xCoords)]
  }
  
  # Note the bg is to be applied to each point
  if(length(bg) != length(xCoords)){
    bg <- rep(bg, ceiling(length(xCoords)/length(bg)))[1:length(xCoords)]
  }
  
  
  # Check that the input data are in the correct format #
  
  
  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addPoints() The vectors containing the X and Y coodinates must be the same length.")
  }
  
  
  # Get the axis limits #
  
  
  # Get the axis limits
  axisLimits <- graphics::par("usr")
  
  
  # Check for NA coordinates #
  
  
  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){
    
    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")
    
    # Check for each of the parameters that can have multiple parameters
    if(length(col.line) == length(xCoords)){
      col.line <- col.line[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty <- lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd <- lwd[-indicesOfNAs]
    }
    if(length(cex) == length(xCoords)){
      cex <- cex[-indicesOfNAs]
    }
    if(length(bg) == length(xCoords)){
      bg <- bg[-indicesOfNAs]
    }
    
    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]
  }
  
  
  # Check if axes are logged #
  
  
  # Check X axis
  xAxisLogged <- FALSE
  if(graphics::par("xlog")){
    
    # Note that X axis was logged
    xAxisLogged <- TRUE
    
    # Log the X coordinates
    xCoords <- log10(xCoords)
    
    # Reset the X axis logged flag - fools points and polygon commands below
    graphics::par(xlog=FALSE)
  }
  
  # Check Y axis
  yAxisLogged <- FALSE
  if(graphics::par("ylog")){
    
    # Note that Y axis was logged
    yAxisLogged <- TRUE
    
    # Log the Y coordinates
    yCoords <- log10(yCoords)
    
    # Reset the Y axis logged flag - fools points and polygon commands below
    graphics::par(ylog=FALSE)
  }
  
  #
  # Store the point information #
  #
  
  # Calculate the height and width of point on current plot
  pointSize <- calculatePointSize(axisLimits, sizeFactor=avoidFactor)
  
  # Store the input coordinates and labels
  # !Note need to make addTextLabels have multiple cex values!
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "N"=length(xCoords), "Heights"=pointSize[1]*cex, 
                    "Widths"=pointSize[2]*cex, "cex"=1)
  

  # Produce a list of alternative locations #

  
  # Generate the alternative locations
  alternativeLocations <- generateAlternativeLocations(axisLimits)
  
  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)
  

  # Create a list to store the information about plotted points #

  
  # Initialise the list to store the information about plotted labels
  plottedPointInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)
  
  
  # Add labels to plot assigning new locations where necessary #
  
  
  # Plot the point label
  for(i in seq_len(pointInfo$N)){
    
    # Set the line characteristics
    lineColour <- setOption(options=col.line, index=i)
    lineType <- setOption(options=lty, index=i)
    lineWidth <- setOption(options=lwd, index=i)
    
    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]
    
    # Get a new location
    newLocationIndex <- chooseNewLocation(pointInfo, i, alternativeLocations, distances,
                                          plottedPointInfo, axisLimits, keepInside)
    
    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 && 
       (tooClose(x, y, height, width, plottedPointInfo) || 
        outsidePlot(x, y, height, width, axisLimits))){
      
      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]
      
      # Add line back to previous location - from the outside of the circle
      graphics::points(x=c(altX, x), y=c(altY, y), type="l", col=col.line, lty=lty, lwd=lwd, xpd=TRUE)
      
      # Add point
      graphics::points(x=altX, y=altY, cex=cex[i], bg=bg[i], ...)
      
      # Append the plotted point information
      plottedPointInfo <- addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedPointInfo)
      
      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]
      
    }else{
      
      # Add point
      graphics::points(x=x, y=y, cex=cex[i], bg=bg[i], ...)
      
      # Append the plotted point information
      plottedPointInfo <- addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedPointInfo)
    }
  }
  
  
  # Return axes logged flags to original state - for if person makes any future plots #
  
  
  graphics::par(xlog=xAxisLogged)
  graphics::par(ylog=yAxisLogged)
  
}

#' Calculate the size of a point on the current plot
#'
#' Function used by \code{addPoints()}
#' @param axisLimits The limits of the X and Y axis: (\code{c(xMin, xMax, yMin, yMax)})
#' @param sizeFactor A number that increases (values > 1) or decreases (values < 1) the amount of space alloted to each point
#' @keywords internal
#' @return Returns a vector containing the width and height of a point
calculatePointSize <- function(axisLimits, sizeFactor=1){
  
  # Get the plotting window size in inches
  plotSizeInches <- graphics::par()$pin # width, height
  widthInches <- plotSizeInches[1]
  heightInches <- plotSizeInches[2]
  
  # Get the plotting window size in the plotting units
  widthX <- axisLimits[2] - axisLimits[1]
  heightY <- axisLimits[4] - axisLimits[3]
  
  # Calculate the size of a point in the current plot
  # Cex=1 is 1/72 inches (https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/pdf.html)
  # Dividing by 72 is far too small - decided to choose 15?!?!
  pointWidth <- (widthX / widthInches) / (15/sizeFactor)
  pointHeight <- (heightY / heightInches) / (15/sizeFactor)
  
  return(c(pointWidth, pointHeight))
}

#' A function to assign value if multiple options are available, can recycle if index is > number of options available
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param colours A single option of vector of options
#' @param index The current index of a label
#' @keywords internal
#' @return Returns the selected option based upon the index provided
setOption <- function(options, index){
  
  # Check if option is null
  option <- NULL
  if(is.null(options) == FALSE){
    # Calculate modulus - the remainder when the index is divided by the number of options provided
    modulus <- index %% length(options)
    
    # Check if modulus non-zero - there is a remainder
    if(modulus != 0){
      
      # Assign option using modulus as index
      option <- options[modulus]
      
      # If no remainder, then index must be the length of the options vector
    }else{
      option <- options[length(options)]
    }
  }
  
  return(option)
}

#' Add the information associated with a text label that has been plotted
#'
#' Function used by \code{addTextLabels()}
#' @param x X coordinate of point of interest
#' @param y Y coodrinate of point of interest
#' @param height The height of the label associated with the point of interest
#' @param width The width of the label associated with the point of interest
#' @param plottedLabelInfo The coordinates and label information about the locations where a label has already plotted
#' @keywords internal
#' @return Returns a list containing information for all the plotted labels, included the one just added
addPlottedLabel <- function(x, y, height, width, plottedLabelInfo){
  
  plottedLabelInfo$X[plottedLabelInfo$N + 1] <- x
  plottedLabelInfo$Y[plottedLabelInfo$N + 1] <- y
  plottedLabelInfo$Heights[plottedLabelInfo$N + 1] <- height
  plottedLabelInfo$Widths[plottedLabelInfo$N + 1] <- width
  
  plottedLabelInfo$N <- plottedLabelInfo$N + 1
  
  return(plottedLabelInfo)
}

#' Plot line from new alternative location back to original
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param altX The X coordinate of new location
#' @param altY The Y coordinate of new location
#' @param x The X coordinate of original location
#' @param y The Y coordinate of original location
#' @param label The label to be plotted. Required to work out when line ends
#' @param cex The number used to scale the size of the label. Required to work out when line ends
#' @param col Colour of line to be plotted
#' @param lty A number detailing the type of line to be plotted. 0: blank, 1: solid, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, and 6: twodash.
#' @param lwd A number to scale the size of plotted line.
#' @param heightPad Multiplyer for label height should added to label to be used to pad height
#' @param widthPad Multiplyer for label width should added to label to be used to pad width
#' @keywords internal
addLineBackToOriginalLocation <- function(altX, altY, x, y, label, cex, col, lty, lwd, heightPad, widthPad){
  
  # Calculate the label width and height
  labelHeight <- graphics::strheight(label, cex=cex)
  labelWidth <- graphics::strwidth(label, cex=cex)
  
  # Calculate amount outer left/right and above/below
  xHalf <- labelWidth * (0.5 + (0.5 * widthPad))
  yHalf <- labelHeight * (0.5 + (0.5 * heightPad))
  
  # Create a set of points marking the boundaries of the label
  xMarkers <- c(seq(from=altX - xHalf, to=altX + xHalf, by=0.05*labelWidth), altX + xHalf)
  yMarkers <- c(seq(from=altY - yHalf, to=altY + yHalf, by=0.05*labelHeight), altY + yHalf)
  
  # Calculate the closest pair of X and Y coordinates to the origin
  closestX <- xMarkers[which.min(abs(xMarkers - x))]
  closestY <- yMarkers[which.min(abs(yMarkers - y))]
  
  # Plot the line
  graphics::points(x=c(closestX, x), y=c(closestY, y), type="l", col=col, lty=lty, lwd=lwd, xpd=TRUE)
}

#' Calculate the heights and widths of the labels in the current plotting window
#'
#' Function used by \code{addTextLabels()}
#' @param pointInfo A list storing the coordinates and labels of input points
#' @param cex The number used to scale the size of the label and therefore its height and width
#' @param heightPad Multiplyer for label height should added to label to be used to pad height
#' @param widthPad Multiplyer for label width should added to label to be used to pad width
#' @keywords internal
#' @return Returns a list storing the coordinates, labels, and the heights and widths of the labels, for input points
calculateLabelHeightsAndWidths <- function(pointInfo, cex, heightPad, widthPad){
  
  # Get the text label heights and lengths
  textHeights <- graphics::strheight(pointInfo$Labels)
  textWidths <- graphics::strwidth(pointInfo$Labels)
  
  # Multiply by cex
  textHeights <- textHeights * cex
  textWidths <- textWidths * cex
  
  # Add padding to widths and heights
  # Note multiplies padding by 2 - stops background polygons being directly adjacent
  pointInfo[["Heights"]] <- textHeights + (2 * heightPad * textHeights)
  pointInfo[["Widths"]] <- textWidths + (2 * widthPad * textWidths)
  
  return(pointInfo)
}

#' Generate a set of alternative locations where labels can be plotted if they overlap with another label
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param axisLimits The limits of the X and Y axis: (\code{c(xMin, xMax, yMin, yMax)})
#' @keywords internal
#' @return Returns a list containing the coordinates of the alternative locations
generateAlternativeLocations <- function(axisLimits){
  
  # Initialise a list to store the alternative locations
  alternativeLocations <- list("X"=c(), "Y"=c())
  
  # Define the spacer for each axis
  spacerX <- 0.01 * (axisLimits[2] - axisLimits[1])
  spacerY <- 0.01 * (axisLimits[4] - axisLimits[3])
  
  # Generate the set of points based upon the spacer
  for(i in seq(axisLimits[1], axisLimits[2], spacerX)){
    for(j in seq(axisLimits[3], axisLimits[4], spacerY)){
      
      alternativeLocations$X[length(alternativeLocations$X) + 1] <- i
      alternativeLocations$Y[length(alternativeLocations$Y) + 1] <- j
    }
  }
  #graphics::points(alternativeLocations$X, alternativeLocations$Y, col=rgb(0,0,0, 0.5), pch=20, xpd=TRUE)
  
  # Note the number of alternative locations created
  alternativeLocations[["N"]] <- length(alternativeLocations$X)
  
  return(alternativeLocations)
}

#' Plot a label with optional polygon background
#'
#' Function used by \code{addTextLabels()}
#' @param x The X coordinate at which label is to be plotted
#' @param y The Y coordinate at which label is to be plotted
#' @param label The label to be plotted
#' @param cex The number used to scale the size of the label
#' @param col The colour of the label to be plotted
#' @param bg The colour of the polygon to be plotted. If NULL no polygon plotted
#' @param border The colour of the polygon border. If NA, no border plotted
#' @param heightPad Multiplyer for label height should added to label to be used to pad height
#' @param widthPad Multiplyer for label width should added to label to be used to pad width
#' @keywords internal
addLabel <- function(x, y, label, cex, col, bg, border, heightPad, widthPad){
  
  # Add a background polygon - if requested
  if(is.null(bg) == FALSE){
    
    # Calculate the height and width of the label
    labelHeight <- graphics::strheight(label, cex=cex)
    labelWidth <- graphics::strwidth(label, cex=cex)
    
    # Calculate amount outer left/right and above/below
    xHalf <- labelWidth * (0.5 + (0.5 * widthPad))
    yHalf <- labelHeight * (0.5 + (0.5 * heightPad))
    
    # Plot the background polygon
    graphics::polygon(x=c(x - xHalf, x - xHalf, x + xHalf, x + xHalf),
                      y=c(y - yHalf, y + yHalf, y + yHalf, y - yHalf),
                      col=bg, border=border, xpd=TRUE)
  }
  
  
  # Add label
  graphics::text(x=x, y=y, labels=label, xpd=TRUE, cex=cex, col=col, adj = c(0.5, 0.75))
}

#' Remove coordinates of alternative locations that are too close to coordinates
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param altXs A vector of X coordinates for alternative locations
#' @param altYs A vector of Y coordinates for alternative locations
#' @param index The index of the point of interest in the coordinate vectors
#' @param textHeight The height of the label to be plotted at the point of interest
#' @param textWidth The width of the label to be plotted at the point of interest
#' @param distances The distances between the actual and alternative locations
#' @keywords internal
#' @return Returns a list of the coordinates of the alternative locations that weren't too close and the distance matrix of the alternate locations to the actual locations
removeLocationAndThoseCloseToItFromAlternatives <- function(altXs, altYs, index, textHeight, textWidth, distances){
  remove <- c(index)
  for(i in 1:length(altXs)){
    
    if(i == index){
      next
    }
    
    if(abs(altXs[index] - altXs[i]) < textWidth &&
       abs(altYs[index] - altYs[i]) < textHeight){
      remove[length(remove) + 1] <- i
    }
  }
  
  altXs <- altXs[-remove]
  altYs <- altYs[-remove]
  distances <- distances[, -remove]
  
  return(list("X" = altXs, "Y" = altYs, "distances"=distances))
}

#' A function to choose (from the alternative locations) a new location for a label to be plotted at
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param pointInfo A list storing the information for the input points
#' @param index The index of the point of interest
#' @param alternativeLocations The coordinates of the alternative locations
#' @param distances The distances between the alternative locations and the input points
#' @param plottedLabelInfo The coordinates and label information about the locations where a label has already plotted
#' @param axisLimits The limits of the X and Y axis: (\code{c(xMin, xMax, yMin, yMax)})
#' @param keepLabelsInside A logical variable indicating whether the labels shouldn't be plotted outside of plotting region
#' @keywords internal
#' @return Returns the index of the chosen alternative location
chooseNewLocation <- function(pointInfo, index, alternativeLocations, distances, plottedLabelInfo, axisLimits, keepLabelsInside){
  
  # graphics::points(alternativeLocations$X, alternativeLocations$Y, pch=19, xpd=TRUE,
  #        col=rgb(1,0,0, distances[index, ] / max(distances[index, ])))
  
  # Get the information about the current point
  x <- pointInfo$X[index]
  y <- pointInfo$Y[index]
  height <- pointInfo$Heights[index] * pointInfo$cex
  width <- pointInfo$Widths[index] * pointInfo$cex
  
  # Get the indices of the alternative locations as an ordered
  orderedAlternateLocationIndices <- order(distances[index, ])
  
  # Initialise a variable to store the index of the selected alternative location
  indexOfSelectedAlternativeLocation <- -1
  
  # Examine each of the alternate locations in order
  for(i in orderedAlternateLocationIndices){
    
    # Get the coordinates of the current alternative location
    altX <- alternativeLocations$X[i]
    altY <- alternativeLocations$Y[i]
    
    # Check current alternative location isn't too close to plotted labels or the plotted input points or label will overlap with plot edges
    if(overlapsWithPlottedPoints(x=altX, y=altY, height=height, width=width, pointInfo=pointInfo) == FALSE &&
       tooClose(x=altX, y=altY, height=height, width=width, plottedLabelInfo) == FALSE &&
       (keepLabelsInside == FALSE || outsidePlot(x=altX, y=altY, height=height, width=width, axisLimits=axisLimits) == FALSE)){
      
      # Store the current index
      indexOfSelectedAlternativeLocation <- i
      break
    }
  }
  
  return(indexOfSelectedAlternativeLocation)
}

#' Checks whether a point is too close to any of the plotted points
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param x X coordinate of point of interest
#' @param y Y coodrinate of point of interest
#' @param height The height of the label associated with the point of interest
#' @param width The width of the label associated with the point of interest
#' @param pointInfo A list storing the information for the input points - that have been plotted
#' @keywords internal
#' @return Returns a logical variable to indicate whether the point of interest was too close to any plotted points
overlapsWithPlottedPoints <- function(x, y, height, width, pointInfo){
  
  result <- FALSE
  for(i in seq_len(pointInfo$N)){
    
    if(abs(x - pointInfo$X[i]) < width && abs(y - pointInfo$Y[i]) < height){
      result <- TRUE
      break
    }
  }
  
  return(result)
}

#' Checks whether adding a label at the current point will ending up being outside of the plotting window
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param x X coordinate of point of interest
#' @param y Y coodrinate of point of interest
#' @param height The height of the label associated with the point of interest
#' @param width The width of the label associated with the point of interest
#' @param axisLimits The limits of the X and Y axis: (\code{c(xMin, xMax, yMin, yMax)})
#' @keywords internal
#' @return Returns a logical variable to indicate whether the point of interest was too close to any plotted labels
outsidePlot <- function(x, y, height, width, axisLimits){
  
  # Calculate half width and height
  halfWidth <- 0.5 * width
  halfHeight <- 0.5* height
  
  # Check if adding a label at the current point would overlap with the plotting window edges
  result <- FALSE
  if(x + halfWidth > axisLimits[2] ||
     x - halfWidth < axisLimits[1] ||
     y + halfHeight > axisLimits[4] ||
     y - halfHeight < axisLimits[3]){
    result <- TRUE
  }
  
  return(result)
}

#' Checks whether a point is too close to any of the plotted labels
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param x X coordinate of point of interest
#' @param y Y coodrinate of point of interest
#' @param height The height of the label associated with the point of interest
#' @param width The width of the label associated with the point of interest
#' @param plottedLabelInfo The coordinates and label information about the locations where a label has already plotted
#' @keywords internal
#' @return Returns a logical variable to indicate whether the point of interest was too close to any plotted labels
tooClose <- function(x, y, height, width, plottedLabelInfo){
  
  # Check if logged axes were used
  if(graphics::par("xlog")){
    x <- log10(x)
  }
  if(graphics::par("ylog")){
    y <- log10(y)
  }
  
  # Check if the current point is too close to any of the plotted locations
  result <- FALSE
  for(i in seq_len(plottedLabelInfo$N)){
    
    if(abs(x - plottedLabelInfo$X[i]) < (0.5 * plottedLabelInfo$Widths[i]) + (0.5 * width) &&
       abs(y - plottedLabelInfo$Y[i]) < (0.5 * plottedLabelInfo$Heights[i]) + (0.5 * height)){
      result <- TRUE
      break
    }
  }
  
  return(result)
}

#' Calculate the euclidean distance between two sets of points. Note: Rescales X axis to match scale of Y
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param pointInfo A list storing the information for the input points
#' @param alternativeLocations A list storing the coordinates of the alternative locations
#' @param axisLimits The limits of the X and Y axis: (\code{c(xMin, xMax, yMin, yMax)})
#' @keywords internal
#' @return Returns the distances between the sets of points provided
euclideanDistancesWithRescaledXAxis <- function(pointInfo, alternativeLocations, axisLimits){
  
  # Calculate the axis ranges
  xRange = axisLimits[2] - axisLimits[1]
  yRange = axisLimits[4] - axisLimits[3]
  
  # Calculate the xFactor
  xFactor <- yRange / xRange
  
  # Initialise a matrix to store distances - note that it is non-symmetric!!!
  distances <- matrix(NA, nrow=pointInfo$N, ncol=alternativeLocations$N)
  
  # Fill the matrix with distances
  for(row in seq_len(nrow(distances))){
    
    for(col in seq_len(ncol(distances))){
      
      # Calculate the distance between the current pair of points
      # REMEMBER to correct the X values for the axes ranges
      distances[row, col] <- euclideanDistance(x1=pointInfo$X[row] * xFactor,
                                               y1=pointInfo$Y[row],
                                               x2=alternativeLocations$X[col] * xFactor,
                                               y2=alternativeLocations$Y[col])
    }
  }
  
  return(distances)
}

#' Calculate the euclidean distance between two points
#'
#' Function used by \code{addTextLabels()} and \code{addPoints()}
#' @param x1 The X coordinate of the first point
#' @param y1 The Y coordinate of the first point
#' @param x2 The X coordinate of the second point
#' @param y2 The Y coordinate of the second point
#' @keywords internal
#' @return Returns the distance between the points provided
euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
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
      abline(h = 0, col = "grey70", lwd = 0.8, lty = 2)
      axis(1, log(c(0.01, 0.1, 1, 10, 100, 1000, 10000)), labels = NA
           , tck = -0.025, col = color.axes, lwd = 0.8, lwd.ticks = 0.8)
      axis(2, trans(yseq), labels = NA
           , tck = -0.025, col = color.axes, lwd = 0.8, lwd.ticks = 0.8)
      
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
           , cex = cex.text, adj = c(0, 0), font = 2, lheight = .6)
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
                         , plocation = "left") {
  
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
    
    # abundance label
    if (length(abundances) > 1) {
      unit = ifelse(abundances[i] > 1, " trees per ha", " tree per ha")
      graphics::text(xlim[2], 0.9*ylim[2]
                     , paste0("Species abundance\n", abundances[i], unit)
                     , adj = 1
                     , cex = 6/7
                     , pos = 2)
    }
    
    # axes
    abline(h = 0, col = "grey70", lty = 2)
    axis(1, cex.axis = 1)
    if (i == 1) axis(2, cex.axis = 1) else axis(2, labels = F)
    
    # inner labels
    if (i %in% labelsx) {
      mtext("absolute latitude (°)", 1, 2.5, cex = 1)
    }
    if (i %in% labelsy) {
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      mtext(paste("stabilizing CNDD", unit), 2, 3.5, las = 0, cex = 1)
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
        if (plocation == "left") {
          loc = 1
          adj = c(0, 1)
        } else {
          loc = which.max(pred_run$x_poly)
          adj = c(1, 1)
        }
        xloc = pred_run$x_poly[loc]
        yloc = pred_run$pred_line[loc] + 0.8*(rev(pred_run$pred_poly)[loc] - pred_run$pred_line[loc])
        yloc = ifelse(yloc < ylim[2], yloc, ylim[2])
        yloc = ifelse(yloc > ylim[1], yloc, ylim[1])
        
        # print p-value
        text(xloc
             , yloc
             , pred_run$pvalue
             , cex = 6/7
             , adj = adj)
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
             , cex = 6/7
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
                          , ylim = NULL) {
  
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
    
    
    # latitude label
    if (is.null(latitude_names)) {
      graphics::text(transAbund(xlim, ref_abund)[2], 0.9*ylim[2], paste0("Latitude\n", latitudes[i], "°")
                     , adj = 1
                     , cex = 6/7
                     , pos = 2)
    } else {
      graphics::text(transAbund(xlim, ref_abund)[2], 0.9*ylim[2], paste0(latitude_names[i], "\n", latitudes[i], "°")
                     , adj = 1
                     , cex = 6/7
                     , pos = 2)
    }
    
    # axes
    abline(h = 0, col = "grey70", lty = 2)
    axis(1, transAbund(c(0.1, 1, 10, 100, 1000, 10000), ref_abund), labels = F)
    axis(1, transAbund(c(0.1, 10, 1000), ref_abund), labels = c(0.1, 10, 1000), tick = F)
    if (latitudes[i] == latitudes[1]) axis(2) else axis(2, labels = F)
    
    
    # inner labels
    if (i %in% labelsx) {
      mtext("abundance (trees per ha)", 1, 2.5, cex = 1)
    }
    if (i %in% labelsy) {
      unit = ifelse(type == "AME", "(% / year)", "(%)")
      mtext(paste("stabilizing CNDD", unit), 2, 3.5, las = 0, cex = 1)
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
             , cex = 6/7
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
    mtext("abundance (trees per ha)", 1, 1, outer = T, cex = 2/3)
  }
  if (labelsy == "outer") {
    unit = ifelse(type == "AME", "(% / year)", "(%)")
    mtext(paste("stabilizing CNDD", unit), 2, 2.5, outer = T, las = 0, cex = 2/3)
  }
  
}


