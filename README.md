# rnmrfind

This package detects individual NMR peaks and broader "regions of interest" (ROI) for 1D NMR data using a combination of sensitivity and resolution enhancement. For details to come in an upcoming publication.

## Installation

The `rnmrfind` package can be installed directly from GitHub using `devtools`:

```
#!R

library(devtools)
install_github('ssokolen/rnmrfind')
```

The package is still under active development and some functionality may change going forward. Further documentation will be added over time.  

## Basic usage

Basic functionality is limited to the `find_roi` and `find_peaks` functions, both of which require spectral data to be provided as a vector of real chemical shift values (assumed to be ppm) as well as a vector of complex (real and imaginary) spectral intensities. The sweep frequency is also needed for internal conversion between Hz and ppm. Example data for ten spectra downloaded from HMDB (https://hmdb.ca/) has been provided in the `nmr.spectra` list.

```
#!R

library(rnmrfind)

# Pull out example fructose data
d <- nmr.spectra$fructose

# Peak detection technically makes use of ROI detection,
# but both can be called separately
roi <- find_roi(d$chemical.shift, d$intensity, sf = 500.13)
print(roi)

peaks <- find_peaks(d$chemical.shift, d$intensity, sf = 500.13)
print(peaks)
```
