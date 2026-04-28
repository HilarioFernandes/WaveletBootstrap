# =============================================================================
# 2_Bootstrap_methods.R
# =============================================================================
# Purpose  : Core bootstrap and wavelet library.
# Chapter  : Multiple
# Inputs   : None
# Outputs  : Functions `block_boot`, `nonboundary_squared_MODWT_coeffs`, `wv_estimates`, `all_squared_MODWT_coeffs`, `wv_estimates_all`, `bootstrap_wavelet`, `distance_cdfs`
# Depends  : waveslim, stats
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running

# calculating wavelet transforms
# install.packages("waveslim")
library(waveslim)

# ecdf()
library(stats)

################################################################################

#' Block Bootstrap Resampling
#' @param N Time series length
#' @param type "MBB" (Moving Block Bootstrap), "NBB" (Non-overlapping Block Bootstrap), or "SB" (Stationary Bootstrap)
#' @param param Expected block length or parameter, depending on `type`. For SB, it is the success probability p. For MBB/NBB, it is the block length.
#' @return Integer vector of length N containing the resampled indexes
block_boot <- function(N, type, param) {
  if (type == "SB") {
    sampled_indexes <- sample(1:N, 1)

    l <- 1
    bern <- 1

    while ((l <= N)) {
      #| (bern == 0)

      bern <- sample(0:1, 1, prob = c(1 - param, param))

      if (bern == 0) {
        sampled_indexes <- c(
          sampled_indexes,
          (sampled_indexes[length(sampled_indexes)] %% N) + 1
        )
      } else {
        sampled_indexes <- c(sampled_indexes, sample(1:N, 1))
      }

      l <- l + 1
    }

    sampled_indexes <- sampled_indexes[1:(length(sampled_indexes) - 1)]
  } else {
    if (type == "MBB") {
      init_indexes <- 1:(N - param + 1)
    }

    if (type == "NBB") {
      init_indexes <- 1 + param * (0:(floor(N / param) - 1))
    }

    amount_blocks <- ceiling(N / param)

    sampled_init_indexes <- sample(init_indexes, amount_blocks, replace = TRUE)

    sampled_indexes <- NULL

    for (index in sampled_init_indexes) {
      sampled_indexes <- c(sampled_indexes, index:(index + param - 1))
    }
  }

  # return(sampled_indexes)
  return(sampled_indexes[1:N])
}

if (FALSE) {
  block_boot(10, "MBB", 2)
  block_boot(10, "NBB", 2)
  block_boot(10, "SB", 0.25)
}
################################################################################

#' Obtain non-boundary squared MODWT wavelet coefficients
#' Useful for obtaining a point estimate of the wavelet variances and bootstrap estimation.
#' Assumes LA8 filter, formula: L_j = (2^j - 1)*(L-1) + 1, with L=8.
#' @param data Numeric vector of time series data
#' @return List of numeric vectors containing squared non-boundary MODWT coefficients per level
nonboundary_squared_MODWT_coeffs <- function(data) {
  N <- length(data)

  # since we will remove boundary coefficients, the total amount of levels is not
  # log2(N), but floor(log2(1+(N-1)/(8-1)))

  n_levels <- floor(log2(1 + (N - 1) / (8 - 1)))

  data.modwt <- modwt(data, n.levels = n_levels)[1:n_levels]

  for (j in 1:n_levels) {
    L_j <- (2^j - 1) * (8 - 1) + 1

    data.modwt[[j]] <- data.modwt[[j]][L_j:length(data.modwt[[j]])]
  }

  coeffs <- lapply(data.modwt, function(x) x^2)

  return(coeffs)
}

if (FALSE) {
  Y <- rnorm(1024)
  Y.coeffs <- nonboundary_squared_MODWT_coeffs(Y)
}

#' Calculate wavelet variance estimate from nonboundary coefficients
#' Assumes LA8 filter
#' @param data Numeric vector of time series data
#' @return Numeric vector of wavelet variance estimates per level
wv_estimates <- function(data) {
  coeffs <- nonboundary_squared_MODWT_coeffs(data)

  estimates <- unlist(lapply(coeffs, function(x) mean(x)))

  return(estimates)
}

if (FALSE) {
  Y.wv <- wv_estimates(Y)
}

#' Obtain all squared MODWT wavelet coefficients
#' @param data Numeric vector of time series data
#' @return List of numeric vectors containing all squared MODWT coefficients per level
all_squared_MODWT_coeffs <- function(data) {
  N <- length(data)

  n_levels <- log2(N)

  data.modwt <- modwt(data, n.levels = n_levels)[1:n_levels]

  coeffs <- lapply(data.modwt, function(x) x^2)

  return(coeffs)
}

if (FALSE) {
  Y.coeffs_all <- all_squared_MODWT_coeffs(Y)
}

#' Calculate the wavelet variance estimate from all coefficients (no boundary adjustment)
#' Assumes LA8 filter
#' @param data Numeric vector of time series data
#' @return Numeric vector of wavelet variance estimates per level
wv_estimates_all <- function(data) {
  coeffs <- all_squared_MODWT_coeffs(data)

  estimates <- unlist(lapply(coeffs, function(x) mean(x)))

  return(estimates)
}

if (FALSE) {
  Y.wv_all <- wv_estimates_all(Y)
}

################################################################################

#' Bootstrap estimates of wavelet variance
#' @param data Numeric vector of time series data
#' @param ordering "bw" (bootstrap then wavelet) or "wb" (wavelet then bootstrap)
#' @param boundlike Logical flag to keep (`TRUE`) or remove (`FALSE`) boundary-like coefficients
#' @param type "MBB", "NBB" or "SB"
#' @param param A function of the time series length N, returning the bootstrap parameter (e.g. `function(N){N^{-1}}`)
#' @param B Number of bootstrap replications
#' @return B x n_levels matrix of bootstrap wavelet variance estimates
bootstrap_wavelet <- function(data, ordering, boundlike, type, param, B) {
  N <- length(data)

  wv_boot_estimates <- NULL

  if (ordering == "bw") {
    n_levels <- floor(log2(1 + (N - 1) / (8 - 1)))

    for (b in 1:B) {
      wv_boot_estimates_temp <- NULL

      indexes_boot <- block_boot(N, type, param(N))

      data_boot <- data[indexes_boot]

      data_boot.coeffs <- nonboundary_squared_MODWT_coeffs(data_boot)[1:n_levels]

      if (boundlike == FALSE) {
        # eliminating boundary-like coefficients
        for (j in 1:length(data_boot.coeffs)) {
          L_j <- (2^j - 1) * (8 - 1) + 1

          for (i in 1:length(data_boot.coeffs[[j]])) {
            if (prod(indexes_boot[i:(i - 1 + L_j)] == indexes_boot[i]:(indexes_boot[i] + L_j - 1)) == 0) {
              data_boot.coeffs[[j]][i] <- NA
            }
          }

          if (all(is.na(data_boot.coeffs[[j]])) == TRUE) {
            wv_boot_estimates_temp <- c(wv_boot_estimates_temp, NA)
          } else {
            wv_boot_estimates_temp <- c(wv_boot_estimates_temp, mean(data_boot.coeffs[[j]], na.rm = TRUE))
          }
        }
      } else {
        for (j in 1:length(data_boot.coeffs)) {
          wv_boot_estimates_temp <- c(wv_boot_estimates_temp, mean(data_boot.coeffs[[j]], na.rm = TRUE))
        }
      }

      wv_boot_estimates <- rbind(wv_boot_estimates, wv_boot_estimates_temp)
    }

    # return(data_boot.coeffs)
  } else if (ordering == "wb") {
    data.coeffs <- nonboundary_squared_MODWT_coeffs(data)

    for (b in 1:B) {
      wv_boot_estimates_temp <- NULL

      for (coeffs in data.coeffs) {
        N_j <- length(coeffs)

        coeffs_boot <- coeffs[block_boot(N_j, type, param(N_j))]

        wv_boot_estimates_temp <- c(wv_boot_estimates_temp, mean(coeffs_boot))
      }

      wv_boot_estimates <- rbind(wv_boot_estimates, wv_boot_estimates_temp)
    }
  }

  rownames(wv_boot_estimates) <- NULL

  return(wv_boot_estimates)
}

if (FALSE) {
  Y <- rnorm(1024)
  wv_boot_estimates <- bootstrap_wavelet(Y, "wb", TRUE, "MBB", function(N) {
    floor(N / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "wb", FALSE, "MBB", function(N) {
    floor(N / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "wb", TRUE, "SB", function(N) {
    N^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "wb", FALSE, "SB", function(N) {
    N^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "bw", TRUE, "MBB", function(N) {
    floor(N / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "bw", FALSE, "MBB", function(N) {
    floor(N / 4)
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "bw", TRUE, "SB", function(N) {
    N^{
      -1
    }
  }, 10)
  wv_boot_estimates <- bootstrap_wavelet(Y, "bw", FALSE, "SB", function(N) {
    N^{
      -1
    }
  }, 10)
}

################################################################################

#' Calculate maximum absolute difference between two empirical CDFs (KS distance)
#' @param data1 First numeric vector
#' @param data2 Second numeric vector
#' @param points Number of points in the approximation grid
#' @return Maximum absolute difference
distance_cdfs <- function(data1, data2, points) {
  m <- min(c(data1, data2), na.rm = TRUE)
  M <- max(c(data1, data2), na.rm = TRUE)

  grid <- seq(m, M, length.out = points)

  ecdf1 <- ecdf(data1)
  ecdf2 <- ecdf(data2)

  diffs <- NULL

  for (x in grid) {
    diffs <- c(diffs, abs(ecdf1(x) - ecdf2(x)))
  }

  return(max(diffs))
}

if (FALSE) {
  data1 <- rnorm(100, 0, 1)
  data2 <- rnorm(100, 1, 1)

  plot(ecdf(data1))
  lines(ecdf(data2))

  distance_cdfs(data1, data2, 1000)
}
