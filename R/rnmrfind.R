#------------------------------------------------------------------------------
#' Find regions of interest (ROI) using sensitivity enhancement
#' 
#' @param chemical.shift Sequence of chemical shift values in ppm
#' @param intensity Relative intensity values corresponding to each chemical
#'                  shift
#' @param sf Sweep frequency
#' @param peak.width Lower end of expected peak widths (in Hz, measured at half
#'                   height)
#' @param dilation.ppm Dilation factor for ROI detection (in ppm)
#' @param noise.threshold Multiple of noise standard deviation to set as cutoff
#'                        between signal and noise.
#' @param R Optional set of decay rates (in Hz) to assess for apodization
#'          (calculated from the data by default)
#' @param format TRUE to output a list of vectors specifying ROI start and stop
#'               (in ppm), FALSE to use indexes (1 to n), NULL to output raw
#'               logic vector of TRUE/FALSE values where TRUE corresponds to an
#'               ROI
#' 
#' @return A vector of TRUE/FALSE values or a list of ROI start, stop vectors in
#'         ppm
#' 
#' @export
find_roi <- function(chemical.shift, intensity, sf, peak.width = 1,
                     dilation.ppm = 0.05, noise.threshold = 4, 
                     R = NULL, format = TRUE) {

  index <- order(chemical.shift)
  chemical.shift <- chemical.shift[index]
  intensity <- intensity[index]

  t <- calculate_time(chemical.shift, sf)
  s <- fft(intensity, inverse = TRUE)

  if ( is.null(R) ) {
    R <- estimate_decay(t, s)
    R <- R  * c(1, 0.75, 0.5)
  }

  n <- estimate_n(chemical.shift, peak.width, 0.5, sf)
  n.dilation <- dilation.ppm/median(diff(chemical.shift))
  rois <- list()

  for ( r in R ) {
    y <- apodize_sensitivity(t, s, r)
    d <- calculate_derivatives(y, n)
    
    pc1 <- combine_pca(d)

    st.dev <- IQR(pc1)/1.35
    roi <- pc1 > ( median(pc1) + noise.threshold * st.dev )

    # Eliminating spurious selections
    roi <- trim(roi, 100)
    roi <- open(roi, max(3, n))
    roi <- dilate(roi, n.dilation)

    rois <- c(rois, list(roi))
  }

  # Select roi selection with most distinct regions
  n.regions <- lapply(rois, function(x) sum(rle(x)$values))
  i <- which.max(n.regions) 
  roi <- rois[[i]]

  if ( is.null(format) ) {
    roi
  } else if ( format ) {
    list_regions(roi, chemical.shift)
  } else {
    list_regions(roi)
  }
}


#------------------------------------------------------------------------------
#' Find peaks using ROI detection and resolution enhancement
#' 
#' @param chemical.shift Sequence of chemical shift values in ppm
#' @param intensity Relative intensity values corresponding to each chemical
#'                  shift
#' @param sf Sweep frequency
#' @param peak.width Lower end of expected peak widths (in Hz, measured at half
#'                   height)
#' @param dilation.ppm Dilation factor for ROI detection (in ppm)
#' @param noise.threshold Multiple of noise standard deviation to set as cutoff
#'                        between signal and noise.
#' @param step.threshold Threshold step size during iterative peak detection
#'                       (smaller numbers increase computational time while
#'                       bigger numbers may miss peaks). Essentially, the noise
#'                       threshold is converted into a percentile value and this
#'                       value is incremented by the step.threshold until it
#'                       reaches 1. Reasonable options for this value are likely
#'                       to fall in the 0.001 to 0.1 range.
#' @param R Optional set of decay rates (in Hz) to assess for apodization
#'          (calculated from the data by default)
#' @param tmax Optional set of tmax parameters to assess for apodization
#' @param format TRUE to output a list of peak locations in ppm, FALSE to output
#'               raw representation of TRUE/FALSE values with TRUE values
#'               specifying peaks
#' 
#' @return A vector of TRUE/FALSE values or a list of peak locations in ppm
#' 
#' @export
find_peaks <- function(chemical.shift, intensity, sf, peak.width = 0.8, 
                       dilation.ppm = 0.01, noise.threshold = 2, 
                       roi.noise.threshold = 4, step.threshold = 0.001, 
                       R = NA, tmax = NA, format = TRUE) {

  index <- order(chemical.shift)
  chemical.shift <- chemical.shift[index]
  intensity <- intensity[index]

  t <- calculate_time(chemical.shift, 500.13)
  s <- fft(intensity, inverse = TRUE)

  if ( is.na(R) ) {
    R <- estimate_decay(t, s)
    R <- R * c(1, 0.75, 0.5)
  }

  param <- expand.grid(R = R, tmax = tmax)

  n <- estimate_n(chemical.shift, peak.width, 0.5, sf)
  peaks <- list()

  # The following requires roi
  roi.mask <- find_roi(chemical.shift, intensity, sf, peak.width, dilation.ppm, 
                       roi.noise.threshold, R, format = NULL)
  roi <- list_regions(roi.mask)
  
  # External loop over apodization parameters
  for ( i in 1:nrow(param) ) {
    r <- param$R[i]
    tmax <- param$tmax[i]

    i <- apodize_resolution(t, s, r, tmax)
    d <- calculate_derivatives(i, n)
    pc1 <- combine_pca(d)

    st.dev <- IQR(pc1)/1.35
    min.cutoff <- median(pc1) + noise.threshold * st.dev

    # Converting cutoff into quantile
    min.cutoff <- sum(pc1 < min.cutoff)/length(pc1)

    # Internal loop over quantile cutoff
    cutoff <- 1

    temp.peaks <- c()

    while ( TRUE ) {
      cutoff <- cutoff - step.threshold
      
      proi <- pc1 > quantile(pc1, cutoff)

      proi <- trim(proi, 100)
      proi <- open(proi, max(3, n))

      # Only keep peaks within roi
      proi <- proi & roi.mask

      new.peaks <- list_peaks(proi, i)
      temp.peaks <- unique(c(temp.peaks, new.peaks))

      # Initial iterations may not find sufficient peaks, skip if that's the case
      if ( length(temp.peaks) < 2 ) {
        next 
      }

      temp.peaks <- merge_peaks(temp.peaks, n)

      if ( cutoff <= min.cutoff ) break
    }

    peaks <- c(peaks, list(temp.peaks))
  }

  # Select peak selection with most peaks
  n.regions <- lapply(peaks, length)

  i <- which.max(n.regions) 
  x <- peaks[[i]]

  if ( is.null(format) ) {
    out <- rep(FALSE, length(chemical.shift))
    out[x] <- TRUE
    out
  } else if ( format ) {
    chemical.shift[x]
  } else {
    x
  }
}


merge_peaks <- function(peaks, d) {

  peaks <- sort(peaks)
  tree <- hclust(dist(peaks), method = "single")
  groups <- split(peaks, (cutree(tree, h = d)))

  unlist(lapply(groups, function(x) round(mean(x))))
}
