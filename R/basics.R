#------------------------------------------------------------------------------
#' Convert chemical shift to time based on sweep frequency
#' 
#' @param chemical.shift Chemical shift in ppm
#' @param sf Sweep frequency in MHz
#' 
#' @return  A vector of time in seconds
#' 
calculate_time <- function(chemical.shift, sf) {

  td <- length(chemical.shift)
  sw <- max(chemical.shift) - min(chemical.shift)
  sw.herz <- sw * sf
  at <- td/sw.herz

  seq(0, at, length.out = td)
}


#------------------------------------------------------------------------------
#' Estimate average decay rate
#' 
#' @param time Acquisition time in seconds.
#' @param signal Relative intensity values in the time domain.
#' @param sf Sweep frequency in MHz
#' 
#' @return Average decay rate for the spectra.
#'
estimate_decay <- function(time, signal, debug = FALSE) {

  # Generate cumulative envelope
  d <- data.frame(x = time, y = Re(signal))
  d$y <- cumsum(abs(d$y))
  d$y <- d$y/max(d$y)

  # Linearize fit to get initial values
  m <- lm(log(1-y) ~ x, data = d[d$y < 0.99, ]) 
  a <- -coef(m)[2] # slope
  b <- exp(coef(m)[1])

  if ( debug ) {
    plot(d$x, d$y, xlab = "Time (s)", ylab = "Cumulative signal")
    lines(d$x, 1 - b * exp(-a * d$x), col = "red")
  }
  
  # Nonlinear least squares to touch up
  m <- nls(y ~ 1 - b * exp(-a * x), data = d[d$y < 0.99, ], 
           start = list(a = a, b = 1))

  a <- coef(m)[1]
  b <- coef(m)[2]
  
  if ( debug ) {
    lines(d$x, 1 - b * exp(-a * d$x), col = "blue", lwd = 2)
    cat(sprintf("Decay rate: %.4f\n", a))
  }
  
  a
}


#------------------------------------------------------------------------------
#' Estimate peak window size for SG filter
#' 
#' @param chemical.shift Chemical shift in ppm
#' @param width Narrowest peak width (at half height) in Hz
#' @param height Fraction of maximum peak height to include
#' @param sf Sweep frequency in MHz
#' 
#' @return Minimum peak widths in points.
#'
estimate_n <- function(chemical.shift, width, height, sf) {

  # Take half of the width for the calculation
  width <- width/2
  
  # Position on x domain where peak takes value of height
  x <- width * sqrt((1 - height)/height)
  
  # Including both sides of peak
  width <- 2 * x
  
  # Converting from width in Hz to width in points
  td <- length(chemical.shift)
  sw <- max(chemical.shift) - min(chemical.shift)
  sw.herz <- sw * sf
  at <- td/sw.herz

  n <- width * at
  
  # Ensuring odd result
  2 * floor(n/2) + 1
}
