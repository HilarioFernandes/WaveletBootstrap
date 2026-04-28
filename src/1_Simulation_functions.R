# =============================================================================
# 1_Simulation_functions.R
# =============================================================================
# Purpose  : Defines time series simulation Models A-D and their ACVS formulas.
# Chapter  : Multiple
# Inputs   : None
# Outputs  : Functions for ACVS (`acvs_*_fun`), data simulation (`Model_*_sim`), and plots.
# Depends  : None
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running

# simulating Gaussian stationary time series
# install.packages("ltsa")
library(ltsa)

# install.packages("fracdiff")
library(fracdiff)

# Set and create output directory for plots
library(grDevices)
OUTPUT_PATH <- file.path(BASE_PATH, "Plots/Plots_modelos")
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)


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
acvs_A_fun <- function(k) {
  e <- exp(1)

  return(2 * sin(pi * e^(-0.1 * abs(k)) / 6))
}


acvs_A <- sapply(0:8192, acvs_A_fun)

#' Simulate from Model A
#' @param N Time series length (must be <= 8192)
#' @return Numeric vector of simulated time series of length N
Model_A_sim <- function(N) {
  U <- DHSimulate(N, acvs_A[1:N])

  Y <- 2 * sqrt(3) * (pnorm(U) - 1 / 2)

  return(Y)
}

if (FALSE) {
  Y <- Model_A_sim(1024)
  plot(Y, type = "l")
}

################################################################################

# Model B

#' Autocovariance sequence for Model B
#' @param k Lag
#' @return Autocovariance at lag k
acvs_B_fun <- function(k) {
  e <- exp(1)

  Sx0 <- log((1 + sqrt(5)) / 2)

  return(log(1 + cos(0.2 * k) * e^(-0.1 * abs(k)) * e^(-Sx0)))
}

acvs_B <- sapply(0:8192, acvs_B_fun)

#' Simulate from Model B
#' @param N Time series length (must be <= 8192)
#' @return Numeric vector of simulated time series of length N
Model_B_sim <- function(N) {
  U <- DHSimulate(N, acvs_B[1:N])

  const <- log((1 + sqrt(5)) / 2)

  Y <- exp(U) - exp(const / 2)

  return(Y)
}

if (FALSE) {
  Y <- Model_B_sim(1024)
  plot(Y, type = "l")
}

################################################################################

# Model C

#' Autocovariance sequence for Model C
#' @param k Lag
#' @return Autocovariance at lag k
acvs_C_fun <- function(k) {
  return(0.9^abs(k))
}

#' Simulate from Model C
#' @param N Time series length
#' @return Numeric vector of simulated time series of length N
Model_C_sim <- function(N) {
  return(arima.sim(
    n = N, list(ar = c(0.9)),
    sd = sqrt(1 - 0.9^2)
  ))
}

if (FALSE) {
  Y <- Model_C_sim(1024)
  plot(Y, type = "l")
}

################################################################################

# Model D

#' Autocovariance sequence for Model D
#' @param d Fractional differencing parameter
#' @param maxlag Maximum lag to calculate
#' @return Numeric vector of length maxlag + 1 containing the ACVS
acvs_D_fun <- function(d, maxlag) {
  x <- numeric(maxlag + 1)
  x[1] <- gamma(1 - 2 * d) / gamma(1 - d)^2
  for (i in 1:maxlag) {
    x[i + 1] <- ((i - 1 + d) / (i - d)) * x[i]
  }
  return(x)
}

acvs_D <- acvs_D_fun(0.25, 32768)

# Note: The DHSimulate-based Model D block below was superseded by the `fracdiff.sim`
# implementation further down. It remains commented for historical context.
#
# Model_D_sim <- function(N){
#
#   U <- DHSimulate(N,acvs_D[1:N])
#
#   return(U)
#
# }
#
# Y <- Model_D_sim(1024)
#
# plot(Y, type = "l")

################################################################################

# Model D

#' Simulate from Model D
#' @param N Time series length
#' @return Numeric vector of simulated time series of length N
Model_D_sim <- function(N) {
  return(fracdiff.sim(n = N, d = 0.25)$series)
}

if (FALSE) {
  Y <- Model_D_sim(1024)
  plot(Y, type = "l")
}

################################################################################

## ---- Example plots ----
## Set to TRUE to enable plot generation on source, or run the block manually.
if (TRUE) {
  set.seed(1)

  YA <- Model_A_sim(128)
  YB <- Model_B_sim(128)
  YC <- Model_C_sim(128)
  YD <- Model_D_sim(128)


  {
    png(
      file = file.path(OUTPUT_PATH, "Plotmodelos.png"),
      width = 1800, height = 900, res = 210
    )

    par(mfrow = c(2, 2))
    par(
      mar = c(3, 2, 2, 1),
      lwd = 1.5
    )
    par(oma = c(1, 0, 0, 0))
    plot(YA,
      type = "l", main = "(a)", xaxt = "n",
      lwd = 1.5
    )
    axis(1, labels = FALSE)
    plot(YB,
      type = "l", main = "(b)", xaxt = "n",
      lwd = 1.5
    )
    axis(1, labels = FALSE)
    plot(YC,
      type = "l", main = "(c)", xlab = "Índice", mgp = c(2, 1, 0),
      lwd = 1.5
    )
    plot(YD,
      type = "l", main = "(d)", xlab = "Índice", mgp = c(2, 1, 0),
      lwd = 1.5
    )

    dev.off()
  }

  # ACVS plot

  acvs_A_plot <- sapply(0:79, acvs_A_fun)
  acvs_B_plot <- sapply(0:79, acvs_B_fun)
  acvs_C_plot <- sapply(0:79, acvs_C_fun)
  acvs_D_plot <- acvs_D_fun(0.25, 79)

  {
    png(
      file = file.path(OUTPUT_PATH, "Plotacvss.png"),
      width = 1800, height = 900, res = 210
    )

    par(mfrow = c(1, 1))
    par(mar = c(4, 4, 2, 1))
    {
      plot(1:80, acvs_B_plot / acvs_B_plot[1],
        type = "l", col = "grey", lty = "dashed", xlab = "Defasagem",
        ylab = "Autocovariância",
        lwd = 1.5
      )
      lines(1:80, acvs_A_plot / acvs_A_plot[1],
        col = "black", lty = "dashed",
        lwd = 1.5
      )
      lines(1:80, acvs_C_plot / acvs_C_plot[1],
        col = "grey",
        lwd = 1.5
      )
      lines(1:80, acvs_D_plot / acvs_D_plot[1],
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
    }

    dev.off()
  }
}
