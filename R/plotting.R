#------------------------------------------------------------------------------
#' Plot NMR spectrum using plotly
#' 
#' @param chemical.shift Sequence of chemical shift values in ppm.
#' @param intensity Complex intensity values corresponding to each chemical
#'                  shift.
#'
plot_spectrum <- function(chemical.shift, intensity, imaginary = FALSE,
                          peaks = c(), regions = list()) {

  if (! require(plotly)) {
    stop('Install the "plotly" package to generate plots')
  }

  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  f <- list(
    family = "Courier New, monospace",
    size = 16,
    color = "#7f7f7f"
  )

  xaxis <- list(
    title = "Chemical shift (ppm)",
    titlefont = f,
    autorange = "reversed"
  )
  yaxis <- list(
    title = "Intensity",
    titlefont = f
  )

  #---------------------------------------

  x <- chemical.shift
  y <- Re(intensity)
  ymax <- max(y)
  ymin <- min(y)

  p <- plot_ly()

  # Regions are generated first
  rect <- function(p, color = "rgba(100,149,237,0.5)") {
    list(
      type = "rect", 
      y0 = ymin,
      y1 = ymax, 
      yref = "y",
      x0 = p[1], 
      x1 = p[2], 
      xref = "x",
      fillcolor = color,
      line = list(color = color)
    )
  }

  if ( length(regions) > 0 ) {
    regions <- lapply(regions, rect)
  } else {
    regions <- list()
  }

  # Then peaks
  line <- function(p, color = "rgba(128,0,0,0.5)") {
    list(
      type = "line", 
      y0 = y[x == p],
      y1 = ymax, 
      x0 = p, 
      x1 = p, 
      line = list(color = color)
    )
  }

  if ( length(peaks) > 0 ) {
    peaks <- lapply(peaks, line)
  } else {
    peaks <- list()
  }

  # And finally everything is plotted
  p <- p %>%
    layout(shapes = c(regions, peaks)) %>%
    add_trace(x = chemical.shift, y = Re(intensity), 
              color = I("black"), name = I("Real"), 
              type = 'scatter', mode = 'lines') %>%
    layout(legend = legend.opts,
           xaxis = xaxis, yaxis = yaxis)

  if ( imaginary ) {
      p <- p %>%     
        add_trace(x = chemical.shift, y = Im(intensity), 
                  color = I("grey"), name = I("Imaginary"), 
                  type = 'scatter', mode = 'lines')
  }

  p
}

#------------------------------------------------------------------------------
#' Plot overlapping NMR spectra using plotly
#' 
#' @param chemical.shift Sequence of chemical shift values in ppm.
#' @param intensity Complex intensity values corresponding to each chemical
#'                  shift.
#'
plot_spectra <- function(chemical.shift, intensities, scale = FALSE) {

  if (! require(plotly)) {
    stop('Install the "plotly" package to generate plots')
  }

  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  f <- list(
    family = "Courier New, monospace",
    size = 16,
    color = "#7f7f7f"
  )

  xaxis <- list(
    title = "Chemical shift (ppm)",
    titlefont = f,
    autorange = "reversed"
  )
  yaxis <- list(
    title = "Intensity",
    titlefont = f
  )

  #---------------------------------------

  p <- plot_ly() %>%
    layout(legend = legend.opts,
           xaxis = xaxis, yaxis = yaxis)

  for ( label in names(intensities) ) {

    y <- Re(intensities[[label]])
    if ( scale ) y <- y/max(y)

    p <- p %>%
      add_trace(x = chemical.shift, y = y, 
                color = I("black"), name = I(label), 
                type = 'scatter', mode = 'lines')
  }

  p
}

