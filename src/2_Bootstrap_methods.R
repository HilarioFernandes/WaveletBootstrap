# =============================================================================
# 2_Bootstrap_methods.R
# =============================================================================
# Purpose  : Core bootstrap and wavelet library.
# Chapter  : Multiple
# Inputs   : None
# Outputs  : Functions `block_boot`, `nonbnd_sq_modwt_coeffs`,
#            `wv_estimates`, `all_sq_modwt_coeffs`, `wv_estimates_all`,
#            `bootstrap_wavelet`, `distance_cdfs`
# Depends  : waveslim, stats
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

# <- SET THIS before running
base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap"

# calculating wavelet transforms
library(waveslim)

library(stats)

################################################################################

#' Block Bootstrap Resampling
#' @param n Time series length
#' @param type "MBB" (Moving Block Bootstrap), "NBB" (Non-overlapping Block
#'   Bootstrap), or "SB" (Stationary Bootstrap)
#' @param param Expected block length or parameter, depending on `type`. For SB,
#'   it is the success probability p. For MBB/NBB, it is the block length.
#' @return Integer vector of length n containing the resampled indexes
block_boot <- function(n, type, param) {
  if (type == "SB") {
    sampled_indexes <- sample(1:n, 1)

    l <- 1
    bern <- 1

    while ((l <= n)) {
      #| (bern == 0)

      bern <- sample(0:1, 1, prob = c(1 - param, param))

      if (bern == 0) {
        sampled_indexes <- c(
          sampled_indexes,
          (sampled_indexes[length(sampled_indexes)] %% n) + 1
        )
      } else {
        sampled_indexes <- c(sampled_indexes, sample(1:n, 1))
      }

      l <- l + 1
    }

    sampled_indexes <- sampled_indexes[1:(length(sampled_indexes) - 1)]
  } else {
    if (type == "MBB") {
      init_indexes <- 1:(n - param + 1)
    }

    if (type == "NBB") {
      init_indexes <- 1 + param * (0:(floor(n / param) - 1))
    }

    amount_blocks <- ceiling(n / param)

    sampled_init_indexes <- sample(init_indexes, amount_blocks, replace = TRUE)

    sampled_indexes <- NULL

    for (index in sampled_init_indexes) {
      sampled_indexes <- c(sampled_indexes, index:(index + param - 1))
    }
  }

  sampled_indexes[1:n]
}

if (FALSE) {
  block_boot(10, "MBB", 2)
  block_boot(10, "NBB", 2)
  block_boot(10, "SB", 0.25)
}
################################################################################

#' Obtain non-boundary squared MODWT wavelet coefficients
#' Useful for obtaining a point estimate of the wavelet variances and bootstrap
#' estimation.
#' Assumes LA8 filter, formula: l_j = (2^j - 1)*(l-1) + 1, with l=8.
#' @param data Numeric vector of time series data
#' @return List of numeric vectors containing squared non-boundary MODWT
#'   coefficients per level
nonbnd_sq_modwt_coeffs <- function(data) {
  n <- length(data)

  # since we will remove boundary coefficients, the total amount of levels
  # is not log2(n), but floor(log2(1+(n-1)/(8-1)))

  n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

  data_modwt <- modwt(data, n.levels = n_levels)[1:n_levels]

  for (j in 1:n_levels) {
    l_j <- (2^j - 1) * (8 - 1) + 1

    data_modwt[[j]] <- data_modwt[[j]][l_j:length(data_modwt[[j]])]
  }

  coeffs <- lapply(data_modwt, function(x) x^2)

  coeffs
}

if (FALSE) {
  y <- rnorm(1024)
  y_coeffs <- nonbnd_sq_modwt_coeffs(y)
}

#' Calculate wavelet variance estimate from nonboundary coefficients
#' Assumes LA8 filter
#' @param data Numeric vector of time series data
#' @return Numeric vector of wavelet variance estimates per level
wv_estimates <- function(data) {
  coeffs <- nonbnd_sq_modwt_coeffs(data)

  estimates <- unlist(lapply(coeffs, function(x) mean(x)))

  estimates
}

if (FALSE) {
  y_wv <- wv_estimates(y)
}

#' Obtain all squared MODWT wavelet coefficients
#' @param data Numeric vector of time series data
#' @return List of numeric vectors containing all squared MODWT coefficients
#'   per level
all_sq_modwt_coeffs <- function(data) {
  n <- length(data)

  n_levels <- log2(n)

  data_modwt <- modwt(data, n.levels = n_levels)[1:n_levels]

  coeffs <- lapply(data_modwt, function(x) x^2)

  coeffs
}

if (FALSE) {
  y_coeffs_all <- all_sq_modwt_coeffs(y)
}

#' Calculate the wavelet variance estimate from all coefficients
#' (no boundary adjustment)
#' Assumes LA8 filter
#' @param data Numeric vector of time series data
#' @return Numeric vector of wavelet variance estimates per level
wv_estimates_all <- function(data) {
  coeffs <- all_sq_modwt_coeffs(data)

  estimates <- unlist(lapply(coeffs, function(x) mean(x)))

  estimates
}

if (FALSE) {
  y_wv_all <- wv_estimates_all(y)
}

################################################################################

#' Bootstrap estimates of wavelet variance
#' @param data Numeric vector of time series data
#' @param ordering "bw" (bootstrap then wavelet) or "wb" (wavelet then
#'   bootstrap)
#' @param boundlike Logical flag to keep (`TRUE`) or remove (`FALSE`)
#'   boundary-like coefficients
#' @param type "MBB", "NBB" or "SB"
#' @param param A function of the time series length n, returning the bootstrap
#'   parameter (e.g. `function(n){n^{-1}}`)
#' @param b_reps Number of bootstrap replications
#' @return b_reps x n_levels matrix of bootstrap wavelet variance estimates
bootstrap_wavelet <- function(data, ordering, boundlike, type, param, b_reps) {
  n <- length(data)

  wv_boot_estimates <- NULL

  if (ordering == "bw") {
    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    for (b in seq_len(b_reps)) {
      wv_boot_estimates_temp <- NULL

      indexes_boot <- block_boot(n, type, param(n))

      data_boot <- data[indexes_boot]

      data_boot_coeffs <-
        nonbnd_sq_modwt_coeffs(data_boot)[1:n_levels]

      if (boundlike == FALSE) {
        # eliminating boundary-like coefficients
        for (j in seq_along(data_boot_coeffs)) {
          l_j <- (2^j - 1) * (8 - 1) + 1

          for (i in seq_along(data_boot_coeffs[[j]])) {
            idx_seq <- indexes_boot[i:(i - 1 + l_j)]
            expected_seq <- indexes_boot[i]:(indexes_boot[i] + l_j - 1)
            if (prod(idx_seq == expected_seq) == 0) {
              data_boot_coeffs[[j]][i] <- NA
            }
          }

          if (all(is.na(data_boot_coeffs[[j]]))) {
            wv_boot_estimates_temp <- c(wv_boot_estimates_temp, NA)
          } else {
            wv_boot_estimates_temp <- c(
              wv_boot_estimates_temp,
              mean(data_boot_coeffs[[j]], na.rm = TRUE)
            )
          }
        }
      } else {
        for (j in seq_along(data_boot_coeffs)) {
          wv_boot_estimates_temp <- c(
            wv_boot_estimates_temp,
            mean(data_boot_coeffs[[j]], na.rm = TRUE)
          )
        }
      }

      wv_boot_estimates <- rbind(wv_boot_estimates, wv_boot_estimates_temp)
    }
  } else if (ordering == "wb") {
    data_coeffs <- nonbnd_sq_modwt_coeffs(data)

    for (b in seq_len(b_reps)) {
      wv_boot_estimates_temp <- NULL

      for (coeffs in data_coeffs) {
        n_j <- length(coeffs)

        coeffs_boot <- coeffs[block_boot(n_j, type, param(n_j))]

        wv_boot_estimates_temp <- c(wv_boot_estimates_temp, mean(coeffs_boot))
      }

      wv_boot_estimates <- rbind(wv_boot_estimates, wv_boot_estimates_temp)
    }
  }

  rownames(wv_boot_estimates) <- NULL

  wv_boot_estimates
}

if (FALSE) {
  y <- rnorm(1024)
  wv_boot_estimates <- bootstrap_wavelet(y, "wb", TRUE, "MBB", function(n) {
    floor(n / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "wb", FALSE, "MBB", function(n) {
    floor(n / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "wb", TRUE, "SB", function(n) {
    n^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "wb", FALSE, "SB", function(n) {
    n^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "bw", TRUE, "MBB", function(n) {
    floor(n / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "bw", FALSE, "MBB", function(n) {
    floor(n / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "bw", TRUE, "SB", function(n) {
    n^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(y, "bw", FALSE, "SB", function(n) {
    n^{
      -1
    }
  }, 10)
}

################################################################################

#' Calculate maximum absolute difference between two empirical CDFs
#' (KS distance)
#' @param data1 First numeric vector
#' @param data2 Second numeric vector
#' @param points Number of points in the approximation grid
#' @return Maximum absolute difference
distance_cdfs <- function(data1, data2, points) {
  min_val <- min(c(data1, data2), na.rm = TRUE)
  max_val <- max(c(data1, data2), na.rm = TRUE)

  grid <- seq(min_val, max_val, length.out = points)

  ecdf1 <- ecdf(data1)
  ecdf2 <- ecdf(data2)

  diffs <- NULL

  for (x in grid) {
    diffs <- c(diffs, abs(ecdf1(x) - ecdf2(x)))
  }

  max(diffs)
}

if (FALSE) {
  data1 <- rnorm(100, 0, 1)
  data2 <- rnorm(100, 1, 1)

  plot(ecdf(data1))
  lines(ecdf(data2))

  distance_cdfs(data1, data2, 1000)
}
