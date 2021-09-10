

#------------------------------------------------------------------------------
#' Perform sensitivity enhancement
#' 
#' @param time Acquisition time in seconds.
#' @param signal Relative intensity values in the time domain.
#' @param R Decay rate paramter in Hz
#' 
#' @return Apodized and fourier transformed spectrum
#' 
apodize_sensitivity <- function(time, signal, R, debug = FALSE) {
  
  f <- exp(-R * time)
  
  if ( debug ) f
  else fft(signal * f)/length(time)
}


#------------------------------------------------------------------------------
#' Perform resolution enhancement
#' 
#' @param time Acquisition time in seconds.
#' @param signal Relative intensity values in the time domain.
#' @param R Decay rate paramter in Hz
#' @param tmax Location of Gaussian bell in seconds (estimated by default)
#' 
#' @return Apodized and fourier transformed spectrum
#' 
apodize_resolution <- function(time, signal, R, tmax = NA, debug = FALSE) {
  
  # Choose t.max as the the point where .99 of the signal is collected
  if ( is.na(tmax) ) {
    accumulated <- cumsum(abs(signal))
    accumulated <- accumulated/max(accumulated)
    t.max <- time[accumulated > 0.99][1]
  }

  f <- exp(R * time) * exp(-(R * time^2)/(2 * t.max))
  
  if ( debug ) f
  else fft(signal * f)/length(time)
}


#------------------------------------------------------------------------------
#' Calculate derivatives using SG
#' 
#' @param intensity Relative intensity values in the frequency domain
#' @param n Windows size (in points) for SG
#' @param smooth TRUE to take the absolute value and smooth the first derivative
#'               of the real data and second derivative of the imaginary
#' 
#' @return A data.frame containing the 1st and 2nd derivatives of the real and
#'         imaginary components (log-scaled)
calculate_derivatives <- function(intensity, n) {

  # Forcing n to a minimum of 5, minimum to avoid interpolation
  if ( n < 5 ) {
    # Disabling warning for now
    #warning("Increasing minimum filter width to avoid interpolation.")
    n <- 5
  }

  r1 <- sgolayfilt(Re(intensity), p = 3, n = n, m = 1)
  r2 <- sgolayfilt(Re(intensity), p = 3, n = n, m = 2)

  # Log filtering has issues with values on either side of zero
  r1 <- abs(r1)
  r2 <- ifelse(r2 < 0, -r2, 0)

  i1 <- sgolayfilt(Im(intensity), p = 3, n = n, m = 1)
  i2 <- sgolayfilt(Im(intensity), p = 3, n = n, m = 2)

  # Log filtering has issues with values on either side of zero
  i1 <- ifelse(i1 < 0, -i1, 0)
  i2 <- abs(i2)

  components <- list(r1, r2, i1, i2)

  d <- do.call(cbind, components)
  colnames(d) <- c("r1", "r2", "i1", "i2")
  d
}


#------------------------------------------------------------------------------
#' Combine data using PCA
#' 
#' @param intensities A data.frame of intensity columns
#' @param scale TRUE to perform log scaling
#' 
#' @return A vector of binary value specifying putative peaks
combine_pca <- function(intensities, scale = FALSE) {

  # Scaling
  f <- function(x) {
    ifelse(x < 0, -log10(-x + 1), log10(x + 1))
  }

  intensities <- as.matrix(intensities) 
  if ( scale ) intensities <- f(intensities)

  # PCA
  pca <- prcomp(intensities, scale = TRUE, rank = 2)
  pc1 <- pca$x[, 1]

  # Ensure consistent signs (signal should end up on the positive end)
  if ( abs(min(pc1)) > abs(max(pc1)) ) pc1 <- -pc1

  pc1
}


#------------------------------------------------------------------------------
#' Wrapper around mmand::dilate to make it more convenient
#' 
#' @param logic A vector of TRUE/FALSE values
#' @param n Dilation width
#' 
#' @return A vector of TRUE/FALSE values
#'
dilate <- function(logic, n) {
  
  k <- c(rep(1, times = (2 * floor(n/2) + 1)))
  as.logical(mmand::dilate(logic, k))
}


#------------------------------------------------------------------------------
#' Wrapper around mmand::opening to make it more convenient
#' 
#' @param logic A vector of TRUE/FALSE values
#' @param n Opening width
#' 
#' @return A vector of TRUE/FALSE values
#'
open <- function(logic, n) {

  k <- c(rep(1, times = (2 * floor(n/2) + 1)))
  as.logical(mmand::opening(logic, k))
}


#------------------------------------------------------------------------------
#' Set n values on either side of vector to FALSE
#' 
#' @param logic A vector of TRUE/FALSE values
#' @param n Number of values to make FALSE on either side of vector
#' 
#' @return A vector of TRUE/FALSE values
#'
trim <- function(logic, n) {

  logic[1:n] <- FALSE
  logic[(length(logic) - n + 1):length(logic)] <- FALSE
  logic
}


#------------------------------------------------------------------------------
#' Convert a series of TRUE/FALSE values into start/stop indexes or values
#' 
#' @param logic A vector of TRUE/FALSE values
#' @param values Optional vector of values to output instead of indexes
#' 
#' @return A list of start/stop pairs (either as indexes or values)
#'
list_regions <- function(logic, values = NULL) {

  n <- length(logic)

  all.indexes <- 1:n
  starts <- all.indexes[c(logic[1], diff(logic) == 1)]
  stops <- all.indexes[c(diff(logic) == -1, logic[n])]

  if ( ! is.null(values) ) {
    starts <- values[starts]
    stops <- values[stops]
  }

  mapply(c, starts, stops, SIMPLIFY = FALSE)
}


#------------------------------------------------------------------------------
#' Convert a series of TRUE/FALSE values into peak positions based on local
#' maxima
#' 
#' @param logic A vector of TRUE/FALSE values
#' @param reference A set of values to check for local maxima
#' @param values Optional vector of values to output instead of indexes
#' 
#' @return A vector of peak positions (either as indexes or values)
#'
list_peaks <- function(logic, reference, values = NULL) {

  # First, get regions
  regions <- list_regions(logic)

  # Then find maximum for each region
  f_maximum <- function(x) {
    index <- x[1]:x[2]
    i <- which.max(Re(reference[index]))
    index[i]
  }

  indexes <- unlist(lapply(regions, f_maximum))

  if ( is.null(values) ) indexes
  else values[indexes]
}


#------------------------------------------------------------------------------
#' Count the number of peaks in each region
#' 
#' @param regions A list of start/stop pairs indicating region bounds
#' @param peaks A vector of peak locations
#' 
#' @return A vector of peak counts per region
#'
count_peaks <- function(regions, peaks) {

  f_count <- function(x) {
    sum( (peaks > x[1]) & (peaks < x[2]) )
  }

  unlist(lapply(regions, f_count))
}


#------------------------------------------------------------------------------
#' Group peaks based on region
#' 
#' @param regions A list of start/stop pairs indicating region bounds
#' @param peaks A vector of peak locations
#' 
#' @return A vector of peak counts per region
#'
group_peaks <- function(regions, peaks) {

  f_group <- function(x) {
    peaks[(peaks > x[1]) & (peaks < x[2])]
  }

  lapply(regions, f_group)
}
