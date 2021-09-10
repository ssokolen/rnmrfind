devtools::load_all("../..")

library(cowplot)
library(grid)
library(gridExtra)
library(tidyverse)

#------------------------------------------------------------------------------
f_data <- function(compound, overwrite = FALSE, suffix = NA, ...) {

  if (! is.na(suffix) ) {
    name <- paste(compound, suffix)
  } else {
    name <- compound
  }

  if ( "data.Rdata" %in% list.files(".") ) {
    load("data.Rdata")
  }

  if (! exists("peak.data") ) {
    peak.data = list()
  }

  if (! overwrite ) {
    if ( name %in% names(peak.data) ) {
      return(peak.data[[name]])
    }
  }

  d <- nmr.spectra[[compound]]

  x <- d$chemical.shift
  y <- d$intensity
  t <- calculate_time(x, 500.13)
  s <- fft(y, inverse = TRUE)
  R <- estimate_decay(t, s)

  y2 <- Re(apodize_resolution(t, s, R))

  peaks <- find_peaks(x, y, 500.13, ...)
  d <- data.frame(x = x, y = Re(y)/max(Re(y)), y2 = Re(y2)/max(Re(y2)))

  out <- list(d, peaks)

  peak.data[[name]] <- out
  save(peak.data, file = "data.Rdata")

  out
}

#------------------------------------------------------------------------------
f_plot <- function(d, lower, upper, offset = 0.005, margin = 30, 
                   comparison = list()) {

  peaks <- d[[2]]
  d <- d[[1]]

  heights <-  d$y[d$x %in% peaks]
  d.peaks <- data.frame(x = peaks, y = heights, method = "rnmrfind")

  methods <- c("rnmrfind")
  for ( method in names(comparison) ) {
    d.new <- comparison[[method]]
    
    d.new$method <- method
    d.peaks <- rbind(d.peaks, d.new)
    methods <- c(methods, method)
  }

  d.peaks$method <- factor(d.peaks$method, levels = methods)

  f_stack <- function(d) {
    d <- d[rev(order(d$method)), ]
    d$y <- seq(1, nrow(d))*offset + d$y
    d
  }

  # Stack different methods
  d.peaks <- d.peaks %>%
    mutate(x = sprintf("%.3f", x)) %>%
    group_by(x) %>%
    do(f_stack(.)) %>%
    mutate(x = as.numeric(x)) %>%
    ungroup()

  d <- filter(d, x > lower, x < upper)
  d.peaks <- filter(d.peaks, x > lower, x < upper)

  p <- ggplot(d) +
    geom_line(aes(x = x, y = y)) +
    geom_line(aes(x = x, y = y2), colour = "grey") +
    scale_x_reverse() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme_bw(18) +
    theme(axis.title = element_blank(),
          plot.margin = margin(l = margin))

  if ( length(comparison) == 0 ) {
    p <- p + geom_point(data = d.peaks, aes(x = x, y = y), shape = 8)
  } else {
    p <- p + 
      geom_point(data = d.peaks, aes(x = x, y = y, shape = method)) +
      scale_shape_manual("Method", values = c(8, 0, 1)) +
      theme(legend.position = "bottom")
  }

  p
}

#------------------------------------------------------------------------------
f_labels <- function(p, xlabel, ylabel, hjust = 0.2) {

  par <- gpar(fontface = "bold", fontsize = 18)

  y.grob <- textGrob(ylabel, gp = par, rot = 90)
  x.grob <- textGrob(xlabel, gp = par, hjust = hjust)

  grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob,
                           padding = unit(2, "line")))

}

#------------------------------------------------------------------------------
f_sim <- function(d, lower, upper, margin = 30) {

  x <- d$chemical.shift
  y <- d$intensity
  t <- calculate_time(x, 500.13)
  s <- fft(y, inverse = TRUE)
  R <- estimate_decay(t, s)

  y2 <- Re(apodize_resolution(t, s, R))

  d <- data.frame(x = x, y = Re(y)/max(Re(y)), y2 = Re(y2)/max(Re(y2)))

  p <- ggplot(filter(d, x > lower, x < upper)) +
    geom_line(aes(x = x, y = y - median(y))) +
    geom_line(aes(x = x, y = y2 - median(y2)), colour = "grey") +
    scale_x_reverse() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme_bw(18) +
    xlab("Chemical shift (ppm)") +
    ylab("Relative intensity") +
    theme(plot.margin = margin(l = margin,r = 10))
}

#------------------------------------------------------------------------------
# Simulated data
load("simulation.Rdata")

d.summary <- d %>%
  group_by(sep, width, snr) %>%
  summarize(percent = sum(peaks & roi)/length(peaks)*100,
            sep = paste(sep, "Hz separation"),
            width = as.character(width))

colours = c("gray70", "gray40", "black")

p1 <- ggplot(d.summary, aes(x = snr, y = percent, shape = width, color = width)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ sep) +
  scale_colour_manual("Peak width (Hz)", values = colours) +
  scale_shape_manual("Peak width (Hz)", values = c(6, 1, 2)) +
  scale_x_log10() +
  xlab("Signal to noise ratio") +
  ylab("Fraction detected") +
  theme_bw(18) +
  theme(legend.position = "bottom")

load("simulation_examples.Rdata")
p2 <- f_sim(filter(d.examples, example == "a"), 10.23, 10.27)
p3 <- f_sim(filter(d.examples, example == "b"),  10.23, 10.27)
p4 <- f_sim(filter(d.examples, example == "c"),  10.23, 10.27)

p.bottom <- plot_grid(p2, p3, p4, labels = c( 'B', 'C', 'D'), 
                   label_size = 18, ncol = 3)

p <- plot_grid(p1, p.bottom, labels = c('A', ''), 
               label_size = 18, ncol = 1, rel_heights = c(1.3, 1))

ggsave("fig4.pdf", width = 18, height = 9)

#---------------------------------------
d <- f_data("folic.acid", peak.width = 0.8, noise.threshold = 2)
p1 <- f_plot(d, 1.959, 2.248, 0.01, 5) +
  annotate("text", x = 2.157, y = 0.13, label = "a") +
  annotate("text", x = 2.147, y = 0.105, label = "b")

d <- f_data("folic.acid", suffix = "wide", peak.width = 1, noise.threshold = 2)
p2 <- f_plot(d, 1.959, 2.248, 0.01, 5) +
  annotate("text", x = 2.157, y = 0.13, label = "a") +
  annotate("text", x = 2.147, y = 0.105, label = "b")

p <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 18, ncol = 1)
p <- f_labels(p, "Chemical shift (ppm)", "Relative intensity")

ggsave("fig1.pdf", plot = p, width = 9, height = 9)

#------------------------------------------------------------------------------

d <- f_data("valine", peak.width = 0.8, noise.threshold = 2)
p1 <- f_plot(d, 2.09, 2.43)
p2 <- f_plot(d, 2.192, 2.329)
p3 <- f_plot(d, 2.350, 2.436, 0.0002)

p <- plot_grid(p1, p2, p3, labels = c('A', 'B', "C"), label_size = 18, ncol = 1)
p <- f_labels(p, "Chemical shift (ppm)", "Relative intensity")

ggsave("fig2.pdf", plot = p, width = 9, height = 9)

#---------------------------------------
d <- f_data("fructose", peak.width = 0.8, noise.threshold = 2)

comparison <- list(
  "rNMR" = read_csv('rnmr_fructose.csv'),
  "speaq" = read_csv('speaq_fructose.csv')
)

p1 <- f_plot(d, 3.98, 4.05, 0.01, 5, comp = comparison) + 
  theme(legend.position = "none")
p2 <- f_plot(d, 3.77, 3.87, 0.01, 5, comp = comparison) + 
  theme(legend.position = "none")
p3 <- f_plot(d, 4.05, 4.08, 0.0005, 5, comp = comparison) + 
  theme(legend.position = "none")
p4 <- f_plot(d, 0.6, 2.3, 6e-5, 5, comp = comparison)

p.top <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), 
                   label_size = 18, ncol = 3)

p <- plot_grid(p.top, p4, labels = c('', 'D'), 
               label_size = 18, ncol = 1)

p <- f_labels(p, "Chemical shift (ppm)", "Relative intensity")

ggsave("fig3.pdf", plot = p, width = 18, height = 9)
