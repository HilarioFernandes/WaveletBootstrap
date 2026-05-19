# =============================================================================
# 10_Comparison_Changepoint_functions.R
# =============================================================================
# Purpose  : Changepoint detection machinery (sequential and non-sequential).
# Chapter  : Chapter 3
# Inputs   : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Outputs  : Block index pairs indicating detected variance changes.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running

source(file.path(base_path, "src", "1_Simulation_functions.R"))
source(file.path(base_path, "src", "2_Bootstrap_methods.R"))

#' Check if ratio falls outside confidence interval
#'
#' @param x Vector containing [ratio, lower_bound, upper_bound]
#' @return 1 if rejected (outside CI), 0 otherwise
interval_check <- function(x) {
  if ((x[1] <= x[2]) | (x[1] >= x[3])) {
    return(1)
  } else {
    return(0)
  }
}

#' Compare two blocks of data for variance change
#'
#' @param x Data vector for block 1
#' @param y Data vector for block 2
#' @param approximation Either "F" (F-test) or "boot_quant" (Bootstrap quantiles)
#' @param alpha_0 Significance level
#' @param jmax Maximum wavelet level to consider
#' @param b Number of bootstrap reps (if approximation == "boot_quant")
#' @return 1 if change detected, 0 otherwise
blocks_comparison <- function(x, y, approximation, alpha_0, jmax, b) {
  x_wv <- wv_estimates(x)
  y_wv <- wv_estimates(y)

  ratios <- y_wv / x_wv

  x_y <- c(x, y)

  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  alpha <- alpha_0 / jmax

  if (approximation == "boot_quant") {
    ratios_boot <- NULL

    for (b in 1:b) {
      boot_indexes <- block_boot(n, "SB", 1 / (4 * log2(n)))

      x_y_boot <- x_y[boot_indexes]
      x_boot_wv <- wv_estimates(x_y_boot[1:n1])
      y_boot_wv <- wv_estimates(x_y_boot[(n1 + 1):(n)])

      ratios_boot <- rbind(ratios_boot, array(y_boot_wv / x_boot_wv))
    }

    quants <- t(apply(ratios_boot, 2, function(x) {
      quantile(x, c(alpha / 2, 1 - alpha / 2))
    })[, 1:jmax])

    rejection <- max(apply(cbind(ratios[1:jmax], quants), 1, interval_check))
  } else if (approximation == "F") {
    js <- 1:jmax

    l_js <- (2^(js) - 1) * (8 - 1) + 1

    n_js <- n - l_js + 1

    eta_js <- sapply(n_js / 2^js, function(x) {
      max(x, 1)
    })

    quants <- t(sapply(eta_js, function(x) {
      qf(c(alpha / 2, 1 - alpha / 2), x, x)
    }))

    rejection <- max(apply(cbind(ratios[1:jmax], quants), 1, interval_check))
  }

  return(rejection)
}


if (FALSE) {
  # mult_factor <- 1
  #
  # n <- 128
  # x <- model_b_sim(2*n)
  # y <- mult_factor*x[(n+1):(2*n)]
  # x <- x[1:n]
  #
  # jmax <- 3
  #
  # alpha_0 <- 0.05
  #
  # plot(c(x,y), type = "l")
  #
  # blocks_comparison(x, y, "F", 0.05, 3, NA)
  # blocks_comparison(x, y, "boot_quant", 0.05, 3, 100)
}

################################################################################

#' Detect changepoints in a time series using block comparisons
#'
#' @param data Time series data vector
#' @param block_size Length of each block
#' @param sequential Logical. If TRUE, compares current block to all previous blocks.
#'        If FALSE, only compares to the immediate previous block.
#' @param approximation Either "F" or "boot_quant"
#' @param alpha_0 Significance level
#' @param jmax Maximum wavelet level to consider
#' @param b Number of bootstrap reps
#' @return Matrix of block index pairs indicating where rejections occurred
changepointdetection <- function(data, block_size, sequential, approximation, alpha_0, jmax, b) {
  amount_blocks <- as.integer(length(data) / block_size)

  changepoint_blocks <- NULL

  if (sequential == TRUE) {
    for (index1 in 2:amount_blocks) {
      for (index2 in (index1 - 1):1) {
        x <- data[(block_size * (index1 - 1) + 1):(index1 * block_size)]
        y <- data[(block_size * (index2 - 1) + 1):(index2 * block_size)]

        alpha <- alpha_0 * (2^(-index2)) / (1 - 2^(1 - index1))

        rejection_temp <- blocks_comparison(y, x, approximation, alpha, jmax, b)

        if (rejection_temp == 1) {
          changepoint_blocks <- rbind(changepoint_blocks, c(index2, index1))

          break
        }
      }
    }
  } else if (sequential == FALSE) {
    for (index1 in 2:amount_blocks) {
      index2 <- index1 - 1

      x <- data[(block_size * (index1 - 1) + 1):(index1 * block_size)]
      y <- data[(block_size * (index2 - 1) + 1):(index2 * block_size)]

      rejection_temp <- blocks_comparison(y, x, approximation, alpha_0, jmax, b)

      if (rejection_temp == 1) {
        changepoint_blocks <- rbind(changepoint_blocks, c(index2, index1))
      }
    }
  }

  return(changepoint_blocks)
}

if (FALSE) {
  # #Example with two blocks of size 128
  #

# mult_factor <- 1
#
# n <- 128
# x <- model_b_sim(2*n)
# y <- mult_factor*x[(n+1):(2*n)]
# x <- x[1:n]
#
# jmax <- 3
#
# alpha_0 <- 0.05
#
# plot(c(x,y), type = "l")

  #
  #
  # blocks_comparison(x, y, "F", 0.05, 3, NA)
  # changepointdetection(c(x,y), 128, FALSE, "F", 0.05, 3, NA)
  # changepointdetection(c(x,y), 128, TRUE, "F", 0.05, 3, NA)
  #
  # blocks_comparison(x, y, "boot_quant", 0.05, 3, 100)
  # changepointdetection(c(x,y), 128, FALSE, "boot_quant", 0.05, 3, 100)
  # changepointdetection(c(x,y), 128, TRUE, "boot_quant", 0.05, 3, 100)
  #
  # #Example with three blocks of size 128
  #

# n <- 128
# x <- model_b_sim(3*n)
# z <- x[(2*n+1):(3*n)]
# y <- x[(n+1):(2*n)]
# x <- x[1:n]
#
# jmax <- 3
#
# alpha_0 <- 0.05
#
# plot(c(x,y), type = "l")

  #
  #
  # changepointdetection(c(x,y,z), 128, FALSE, "F", 0.05, 3, NA)
  # changepointdetection(c(x,y,z), 128, TRUE, "F", 0.05, 3, NA)
  #
  # changepointdetection(c(x,y,z), 128, FALSE, "boot_quant", 0.05, 3, 100)
  # changepointdetection(c(x,y,z), 128, TRUE, "boot_quant", 0.05, 3, 100)
}
