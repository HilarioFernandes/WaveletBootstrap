# =============================================================================
# 1_Simulation_functions.R
# =============================================================================
# Purpose  : Defines time series simulation Models A-D and their ACVS formulas.
# Chapter  : Multiple
# Inputs   : None
# Outputs  : Functions for ACVS (`acvs_*_fun`), data simulation (`Model_*_sim`),
#            and plots.
# Depends  : None
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

# <- SET THIS before running
base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap"

# simulating Gaussian stationary time series

library(ltsa)


library(fracdiff)

# Set and create output directory for plots
library(grDevices)
output_path <- file.path(base_path, "Plots/Plots_modelos")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


################################################################################

#' Number Formatting function
#' @param val Numeric value
#' @return Formatted string with exactly two decimal places and no leading zero
numformat <- function(val) {
  sub("^(-?)0.", "\\1.", sprintf("%.2f", val))
}

################################################################################

# Model A

#' Autocovariance sequence for Model A
#' @param k Lag
#' @return Autocovariance at lag k
acvs_a_fun <- function(k) {
  e <- exp(1)

  2 * sin(pi * e^(-0.1 * abs(k)) / 6)
}


acvs_a <- sapply(0:8192, acvs_a_fun)

#' Simulate from Model A
#' @param n Time series length (must be <= 8192)
#' @return Numeric vector of simulated time series of length n
model_a_sim <- function(n) {
  u <- DHSimulate(n, acvs_a[1:n])

  y <- 2 * sqrt(3) * (pnorm(u) - 1 / 2)

  y
}

if (FALSE) {
  y <- model_a_sim(1024)
  plot(y, type = "l")
}

################################################################################

# Model b

#' Autocovariance sequence for Model b
#' @param k Lag
#' @return Autocovariance at lag k
acvs_b_fun <- function(k) {
  e <- exp(1)

  sx0 <- log((1 + sqrt(5)) / 2)

  log(1 + cos(0.2 * k) * e^(-0.1 * abs(k)) * e^(-sx0))
}

acvs_b <- sapply(0:8192, acvs_b_fun)

#' Simulate from Model b
#' @param n Time series length (must be <= 8192)
#' @return Numeric vector of simulated time series of length n
model_b_sim <- function(n) {
  u <- DHSimulate(n, acvs_b[1:n])

  const <- log((1 + sqrt(5)) / 2)

  y <- exp(u) - exp(const / 2)

  y
}

if (FALSE) {
  y <- model_b_sim(1024)
  plot(y, type = "l")
}

################################################################################

# Model C

#' Autocovariance sequence for Model C
#' @param k Lag
#' @return Autocovariance at lag k
acvs_c_fun <- function(k) {
  0.9^abs(k)
}

#' Simulate from Model C
#' @param n Time series length
#' @return Numeric vector of simulated time series of length n
model_c_sim <- function(n) {
  arima.sim(
    n = n, list(ar = c(0.9)),
    sd = sqrt(1 - 0.9^2)
  )
}

if (FALSE) {
  y <- model_c_sim(1024)
  plot(y, type = "l")
}

################################################################################

# Model D

#' Autocovariance sequence for Model D
#' @param d Fractional differencing parameter
#' @param maxlag Maximum lag to calculate
#' @return Numeric vector of length maxlag + 1 containing the ACVS
acvs_d_fun <- function(d, maxlag) {
  x <- numeric(maxlag + 1)
  x[1] <- gamma(1 - 2 * d) / gamma(1 - d)^2
  for (i in 1:maxlag) {
    x[i + 1] <- ((i - 1 + d) / (i - d)) * x[i]
  }
  x
}

acvs_d <- acvs_d_fun(0.25, 32768)



#

#

#

#

#

#


################################################################################

# Model D

#' Simulate from Model D
#' @param n Time series length
#' @return Numeric vector of simulated time series of length n
model_d_sim <- function(n) {
  fracdiff.sim(n = n, d = 0.25)$series
}

if (FALSE) {
  y <- model_d_sim(1024)
  plot(y, type = "l")
}

################################################################################

## ---- Example plots ----
## Set to TRUE to enable plot generation on source, or run the block manually.
set.seed(1)

ya <- model_a_sim(128)
yb <- model_b_sim(128)
yc <- model_c_sim(128)
yd <- model_d_sim(128)


png(
  file = file.path(output_path, "Plotmodelos.png"),
  width = 1800, height = 900, res = 210
)

par(mfrow = c(2, 2))
par(
  mar = c(3, 2, 2, 1),
  lwd = 1.5
)
par(oma = c(1, 0, 0, 0))
plot(ya,
  type = "l", main = "(a)", xaxt = "n",
  lwd = 1.5
)
axis(1, labels = FALSE)
plot(yb,
  type = "l", main = "(b)", xaxt = "n",
  lwd = 1.5
)
axis(1, labels = FALSE)
plot(yc,
  type = "l", main = "(c)", xlab = "Índice", mgp = c(2, 1, 0),
  lwd = 1.5
)
plot(yd,
  type = "l", main = "(d)", xlab = "Índice", mgp = c(2, 1, 0),
  lwd = 1.5
)

dev.off()

# ACVS plot

acvs_a_plot <- sapply(0:79, acvs_a_fun)
acvs_b_plot <- sapply(0:79, acvs_b_fun)
acvs_c_plot <- sapply(0:79, acvs_c_fun)
acvs_d_plot <- acvs_d_fun(0.25, 79)

png(
  file = file.path(output_path, "Plotacvss.png"),
  width = 1800, height = 900, res = 210
)

par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(1:80, acvs_b_plot / acvs_b_plot[1],
  type = "l", col = "grey", lty = "dashed", xlab = "Defasagem",
  ylab = "Autocovariância",
  lwd = 1.5
)
lines(1:80, acvs_a_plot / acvs_a_plot[1],
  col = "black", lty = "dashed",
  lwd = 1.5
)
lines(1:80, acvs_c_plot / acvs_c_plot[1],
  col = "grey",
  lwd = 1.5
)
lines(1:80, acvs_d_plot / acvs_d_plot[1],
  col = "black",
  lwd = 1.5
)
abline(
  h = 0, lty = "dotdash", col = "lightgrey",
  lwd = 1.5
)
legend(70, 1,
  legend = c("(a)", "(b)", "(c)", "(d)"),
  col = c("black", "grey", "grey", "black"),
  lty = c("dashed", "dashed", "solid", "solid")
)

dev.off()
