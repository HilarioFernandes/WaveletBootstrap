# =============================================================================
# 6_char_scales.R
# =============================================================================
# Purpose  : Characteristic scales study (scales where wavelet variance is maximised).
# Chapter  : Chapter 3, Appendix b.3
# Inputs   : None
# Outputs  : Characteristic scale estimates and CIs.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "1_Simulation_functions.R"))
source(file.path(base_path, "src", "2_Bootstrap_methods.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

# Set and create output directory for plots
output_path <- file.path(base_path, "Plots/Plots_6")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

################################################################################

# characteristic scales

#' Find level with maximum wavelet variance
#'
#' @param wv_est Vector of wavelet variance estimates
#' @return Index of the maximum value
max_char_scale <- function(wv_est) {
  return(unname(which.max(wv_est)))
}

if (FALSE) {
  # n <- 2048
  # y <- model_a_sim(n)
  # wv_est <- wv_estimates(y)
  # plot(log2(wv_est))
  #
  # max_char_scale(wv_est)
}


#' Point estimate of largest local characteristic scale (level)
#'
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return Quadratic interpolation of the characteristic scale
char_scale_est <- function(wv_est, j_0) {
  log_wv_est <- log2(wv_est)

  if (j_0 == 1 || j_0 == length(wv_est)) {
    return(-1)
  }

  # we assume j_0 > 1

  x <- (j_0 - 1):(j_0 + 1)

  y <- log_wv_est[(j_0 - 1):(j_0 + 1)]

  quadratic <- solve(cbind(1, x, x^2), y)

  return(-quadratic[2] / (2 * quadratic[3]))
}

if (FALSE) {
  # n <- 2048
  # y <- model_a_sim(n)
  # wv_est <- wv_estimates(y)
  # plot(log2(wv_est))
  #
  # char_scale_est(wv_est,max_char_scale(wv_est))
}


#' beta_hat calculation for quadratic interpolation
#'
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return Vector of beta coefficients
beta_hat_est <- function(wv_est, j_0) {
  H <- rbind(c(-1 / 2, 0, 1 / 2), c(1, -2, 1))

  return(as.vector(H %*% log2(wv_est[(j_0 - 1):(j_0 + 1)])))
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_a_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # wv_est <- wv_estimates(y)
  #
  # beta_hat_est(wv_est, 7)
}

#' Function that calculates the terms \hat{s}
#'
#' @param W List with all wavelet coefficients
#' @param n Time series length
#' @param l Array with the filter lengths
#' @return A list with biased autocovariance estimates for each scale
biased_est_acvs <- function(W, n, l) {
  n_levels <- length(l)

  output <- vector(mode = "list", length = n_levels)

  for (j in 1:n_levels) {
    temp <- NULL

    for (tau in 0:(n - l[j])) {
      sum <- 0

      for (t in (l[j] - 1):(n - tau - 1)) {
        sum <- sum + W[[j]][t] * W[[j]][t + tau]
      }

      temp <- c(temp, sum / (n - l[j] + 1))
    }

    output[[j]] <- temp
  }

  return(output)
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_a_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # coeffs <- modwt(y, n.levels = n_levels)[1:n_levels]
  #
  # biased_est_acvs(coeffs, n, l)
}

#' Function that calculates the matrix \sigma_1 from eq (8)
#'
#' @param s_est List of biased acvs estimates
#' @param wv_est Vector of wavelet variance estimates
#' @param n Time series length
#' @param l Array of filter lengths
#' @param j_0 Target level index
#' @return The sigma1 matrix estimate
sigma1_est <- function(s_est, wv_est, n, l, j_0) {
  temp <- matrix(NA, nrow = 3, ncol = 3)

  for (k_1 in 1:3) {
    for (k_2 in k_1:3) {
      sum <- 0

      for (tau in 1:(n - l[j_0 + (k_2 - 2)])) {
        # the lag window chosen is the same as in the paper

        sum <- sum + s_est[[j_0 + (k_1 - 2)]][tau] * s_est[[j_0 + (k_2 - 2)]][tau]
      }

      sum <- (2 * sum + wv_est[j_0 + (k_1 - 2)] * wv_est[j_0 + (k_2 - 2)]) / 2

      temp[k_2, k_1] <- 2 * sum / (n - l[k_1] + 1)
    }
  }

  return(temp)
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_c_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # j_0 <- 3
  #
  # coeffs <- modwt(y, n.levels = n_levels)[1:n_levels]
  #
  # s_est <- biased_est_acvs(coeffs, n, l)
  #
  # wv_est <- wv_estimates(y)
  #
  # sigma1_est(s_est, wv_est, n, l, j_0)
}


#' Function that calculates \sigma_2 from eq (9)
#'
#' @param sigma1_est The sigma1 matrix estimate
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return The sigma2 matrix estimate
sigma2_est <- function(sigma1_est, wv_est, j_0) {
  temp <- matrix(NA, nrow = 3, ncol = 3)

  for (k_1 in 1:3) {
    for (k_2 in k_1:3) {
      sum <- sigma1_est[k_2, k_1] / (wv_est[k_1 + j_0 - 2] * wv_est[k_2 + j_0 - 2] * log(2)^2)

      sum <- sum + 2 * (sigma1_est[k_1, k_1] * sigma1_est[k_2, k_2] + sigma1_est[k_2, k_1]^2) / (wv_est[k_1 + j_0 - 2]^2 * wv_est[k_2 + j_0 - 2]^2 * log(2)^2)


      temp[k_1, k_2] <- sum

      if (k_1 != k_2) {
        temp[k_2, k_1] <- sum
      }
    }
  }

  return(temp)
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_c_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # j_0 <- 3
  #
  # coeffs <- modwt(y, n.levels = n_levels)[1:n_levels]
  #
  # s_est <- biased_est_acvs(coeffs, n, l)
  #
  # wv_est <- wv_estimates(y)
  #
  # sigma1_est <- sigma1_est(s_est, wv_est, n, l, j_0)
  #
  # sigma2_est(sigma1_est, wv_est, j_0)
}


#' Function that calculates the var-cov matrix associated to \hat{\beta}
#'
#' @param sigma2_est The sigma2 matrix estimate
#' @return The variance-covariance matrix for beta coefficients
var_beta_est <- function(sigma2_est) {
  H <- rbind(c(-1 / 2, 0, 1 / 2), c(1, -2, 1))

  return(H %*% sigma2_est %*% t(H))
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_c_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # j_0 <- 3
  #
  # coeffs <- modwt(y, n.levels = n_levels)[1:n_levels]
  #
  # s_est <- biased_est_acvs(coeffs, n, l)
  #
  # wv_est <- wv_estimates(y)
  #
  # sigma1_est <- sigma1_est(s_est, wv_est, n, l, j_0)
  #
  # sigma2_est <- sigma2_est(sigma1_est, wv_est, j_0)
  #
  # var_beta_est(sigma2_est)
}


#' Function that calculates the variance of \hat{\kappa} from eq (10)
#'
#' @param var_beta_est Variance-covariance matrix of beta
#' @param beta_hat Vector of beta coefficient estimates
#' @return Variance of the characteristic scale estimate
var_kappa <- function(var_beta_est, beta_hat) {
  sum <- var_beta_est[1, 1] / beta_hat[2]^2 + (beta_hat[1]^2 * var_beta_est[2, 2]) / beta_hat[2]^4

  sum <- sum + (var_beta_est[1, 1] * var_beta_est[2, 2] + 2 * var_beta_est[1, 2]) / beta_hat[2]^4

  sum <- sum + (3 * beta_hat[1]^2 * var_beta_est[2, 2]^2) / beta_hat[2]^6

  sum <- sum - (2 * beta_hat[1] * var_beta_est[1, 2]) / beta_hat[2]^3

  return(sum)
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_a_sim(n)
  #
  # n_levels <- floor(log2(1+(n-1)/(8-1)))
  #
  # l <- (2^(1:n_levels)-1)*(8-1)+1
  #
  # coeffs <- modwt(y, n.levels = n_levels)[1:n_levels]
  #
  # s_est <- biased_est_acvs(coeffs, n, l)
  #
  # wv_est <- wv_estimates(y)
  #
  # j_0 <- max_char_scale(wv_est)
  #
  # plot(log2(wv_est))
  #
  # sigma1_est <- sigma1_est(s_est, wv_est, n, l, j_0)
  #
  # sigma2_est <- sigma2_est(sigma1_est, wv_est, j_0)
  #
  # var_beta_est <- var_beta_est(sigma2_est)
  #
  # beta_hat <- beta_hat_est(wv_est, max_char_scale(wv_est))
  #
  # var_kappa(var_beta_est, beta_hat)
}

#' Calculate confidence intervals for characteristic scales
#'
#' @param data Time series data vector
#' @param alpha Significance level
#' @return A vector with [lower, upper] limits
ci_calc <- function(data, alpha) {
  n <- length(data)

  n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

  l <- (2^(1:n_levels) - 1) * (8 - 1) + 1

  coeffs <- modwt(data, n.levels = n_levels)[1:n_levels]

  s_est <- biased_est_acvs(coeffs, n, l)

  wv_est <- wv_estimates(data)

  j_0 <- max_char_scale(wv_est)

  if (j_0 == 1 || j_0 == n_levels) {
    return(c(-1, -1))
  }

  sigma1_est <- sigma1_est(s_est, wv_est, n, l, j_0)

  sigma2_est <- sigma2_est(sigma1_est, wv_est, j_0)

  var_beta_est <- var_beta_est(sigma2_est)

  beta_hat <- beta_hat_est(wv_est, j_0)

  var_kappa <- var_kappa(var_beta_est, beta_hat)

  char_scale_est <- char_scale_est(wv_est, j_0)

  if (is.nan(sqrt(var_kappa))) {
    return(c(-1, -1))
  }

  return(c(
    char_scale_est + qnorm(alpha / 2) * sqrt(var_kappa),
    char_scale_est + qnorm(1 - alpha / 2) * sqrt(var_kappa)
  ))
}

if (FALSE) {
  # n <- 2048
  #
  # y <- model_a_sim(n)
  #
  # ci_calc(y, 0.05)
}

################################################################################

# study for CI

set.seed(45)

b <- if (test_mode) 5 else 100
iterations <- if (test_mode) 2 else 100

char_scale_est_a <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_b <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_c <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_d <- list("128" = NULL, "512" = NULL, "2048" = NULL)

{
  a_wv_cs_ci <- list(
    "128" = c(NULL, NULL),
    "512" = c(NULL, NULL),
    "2048" = c(NULL, NULL)
  )
  b_wv_cs_ci <- list(
    "128" = c(NULL, NULL),
    "512" = c(NULL, NULL),
    "2048" = c(NULL, NULL)
  )
  c_wv_cs_ci <- list(
    "128" = c(NULL, NULL),
    "512" = c(NULL, NULL),
    "2048" = c(NULL, NULL)
  )
  d_wv_cs_ci <- list(
    "128" = c(NULL, NULL),
    "512" = c(NULL, NULL),
    "2048" = c(NULL, NULL)
  )
}

{
  a_wv_sb_2 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  a_wv_sb_4 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  a_wv_sb_8 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  b_wv_sb_2 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  b_wv_sb_4 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  b_wv_sb_8 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  c_wv_sb_2 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  c_wv_sb_4 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  c_wv_sb_8 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  d_wv_sb_2 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  d_wv_sb_4 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  d_wv_sb_8 <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )
}


for (iter in 1:iterations) {
  print(iter)

  for (i in 1:3) {
    # print(c(iter,i))

    n <- 2^((2 * i - 1) + 6)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    # simulating from each model
    ya <- model_a_sim(n)
    yb <- model_b_sim(n)
    yc <- model_c_sim(n)
    yd <- model_d_sim(n)

    # calculating the point estimates of the characteristic scales
    wv_est_a <- wv_estimates(ya)
    wv_est_b <- wv_estimates(yb)
    wv_est_c <- wv_estimates(yc)
    wv_est_d <- wv_estimates(yd)

    char_scale_est_a[[i]] <- c(char_scale_est_a[[i]], char_scale_est(wv_est_a, max_char_scale(wv_est_a)))
    char_scale_est_b[[i]] <- c(char_scale_est_b[[i]], char_scale_est(wv_est_b, max_char_scale(wv_est_b)))
    char_scale_est_c[[i]] <- c(char_scale_est_c[[i]], char_scale_est(wv_est_c, max_char_scale(wv_est_c)))
    char_scale_est_d[[i]] <- c(char_scale_est_d[[i]], char_scale_est(wv_est_d, max_char_scale(wv_est_d)))

    # calculating the confidence intervals via the normal approximation
    a_wv_cs_ci[[i]] <- rbind(a_wv_cs_ci[[i]], ci_calc(ya, 0.05))
    b_wv_cs_ci[[i]] <- rbind(b_wv_cs_ci[[i]], ci_calc(yb, 0.05))
    c_wv_cs_ci[[i]] <- rbind(c_wv_cs_ci[[i]], ci_calc(yc, 0.05))
    d_wv_cs_ci[[i]] <- rbind(d_wv_cs_ci[[i]], ci_calc(yd, 0.05))

    # bootstrap wavelet estimates

    a_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    a_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    a_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    b_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    b_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    b_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    c_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    c_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    c_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    d_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    d_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    d_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)
  }
}

################################################################################

true_char_scales <- function(list1) {
  output <- NULL

  for (i in 1:3) {
    if (length(which(list1[[i]] == -1)) > iterations / 2) {
      output <- c(output, -1)
    } else {
      temp <- list1[[i]]

      output <- c(output, median(temp[temp > 0]))
    }
  }

  return(output)
}

char_scale_a <- true_char_scales(char_scale_est_a)
char_scale_b <- true_char_scales(char_scale_est_b)
char_scale_c <- true_char_scales(char_scale_est_c)
char_scale_d <- true_char_scales(char_scale_est_d)


char_scale_est_boot <- function(list1) {
  temp <- list(
    "128" = vector("list", length = iterations),
    "512" = vector("list", length = iterations),
    "2048" = vector("list", length = iterations)
  )

  temp_fun <- function(x) {
    char_scale_est(x, max_char_scale(x))
  }

  for (i in 1:3) {
    for (iter in 1:iterations) {
      temp[[i]][[iter]] <- apply(list1[[i]][[iter]], 1, temp_fun)
    }
  }

  return(temp)
}

char_scale_est_a_sb_2_boot <- char_scale_est_boot(a_wv_sb_2)
char_scale_est_a_sb_4_boot <- char_scale_est_boot(a_wv_sb_4)
char_scale_est_a_sb_8_boot <- char_scale_est_boot(a_wv_sb_8)

char_scale_est_b_sb_2_boot <- char_scale_est_boot(b_wv_sb_2)
char_scale_est_b_sb_4_boot <- char_scale_est_boot(b_wv_sb_4)
char_scale_est_b_sb_8_boot <- char_scale_est_boot(b_wv_sb_8)

char_scale_est_c_sb_2_boot <- char_scale_est_boot(c_wv_sb_2)
char_scale_est_c_sb_4_boot <- char_scale_est_boot(c_wv_sb_4)
char_scale_est_c_sb_8_boot <- char_scale_est_boot(c_wv_sb_8)

char_scale_est_d_sb_2_boot <- char_scale_est_boot(d_wv_sb_2)
char_scale_est_d_sb_4_boot <- char_scale_est_boot(d_wv_sb_4)
char_scale_est_d_sb_8_boot <- char_scale_est_boot(d_wv_sb_8)

ci_calc_boot <- function(list1, alpha) {
  output <- list(c(NULL, NULL), c(NULL, NULL), c(NULL, NULL))

  for (i in 1:3) {
    for (iter in 1:iterations) {
      if (length(which(list1[[i]][[iter]] == -1)) > b / 2) {
        output[[i]] <- rbind(output[[i]], c(-1, -1))
      } else {
        temp <- list1[[i]][[iter]]

        output[[i]] <- rbind(output[[i]], quantile(temp[temp > 0], c(alpha / 2, 1 - alpha / 2)))
      }
    }
  }

  return(output)
}


a_sb_2_ci_boot <- ci_calc_boot(char_scale_est_a_sb_2_boot, 0.05)
a_sb_4_ci_boot <- ci_calc_boot(char_scale_est_a_sb_4_boot, 0.05)
a_sb_8_ci_boot <- ci_calc_boot(char_scale_est_a_sb_8_boot, 0.05)

b_sb_2_ci_boot <- ci_calc_boot(char_scale_est_b_sb_2_boot, 0.05)
b_sb_4_ci_boot <- ci_calc_boot(char_scale_est_b_sb_4_boot, 0.05)
b_sb_8_ci_boot <- ci_calc_boot(char_scale_est_b_sb_8_boot, 0.05)

c_sb_2_ci_boot <- ci_calc_boot(char_scale_est_c_sb_2_boot, 0.05)
c_sb_4_ci_boot <- ci_calc_boot(char_scale_est_c_sb_4_boot, 0.05)
c_sb_8_ci_boot <- ci_calc_boot(char_scale_est_c_sb_8_boot, 0.05)

d_sb_2_ci_boot <- ci_calc_boot(char_scale_est_d_sb_2_boot, 0.05)
d_sb_4_ci_boot <- ci_calc_boot(char_scale_est_d_sb_4_boot, 0.05)
d_sb_8_ci_boot <- ci_calc_boot(char_scale_est_d_sb_8_boot, 0.05)


ci_calc_boot_perc <- function(list1, array2) {
  output <- NULL

  for (i in 1:3) {
    temp_fun <- function(x) {
      (x[1] <= array2[i]) & (x[2] >= array2[i])
    }

    output <- c(output, sum(apply(list1[[i]], 1, temp_fun)) / iterations)
  }

  return(output)
}

a_sb_2_ci_boot_perc <- ci_calc_boot_perc(a_sb_2_ci_boot, char_scale_a)
a_sb_4_ci_boot_perc <- ci_calc_boot_perc(a_sb_4_ci_boot, char_scale_a)
a_sb_8_ci_boot_perc <- ci_calc_boot_perc(a_sb_8_ci_boot, char_scale_a)

b_sb_2_ci_boot_perc <- ci_calc_boot_perc(b_sb_2_ci_boot, char_scale_b)
b_sb_4_ci_boot_perc <- ci_calc_boot_perc(b_sb_4_ci_boot, char_scale_b)
b_sb_8_ci_boot_perc <- ci_calc_boot_perc(b_sb_8_ci_boot, char_scale_b)

c_sb_2_ci_boot_perc <- ci_calc_boot_perc(c_sb_2_ci_boot, char_scale_c)
c_sb_4_ci_boot_perc <- ci_calc_boot_perc(c_sb_4_ci_boot, char_scale_c)
c_sb_8_ci_boot_perc <- ci_calc_boot_perc(c_sb_8_ci_boot, char_scale_c)

d_sb_2_ci_boot_perc <- ci_calc_boot_perc(d_sb_2_ci_boot, char_scale_d)
d_sb_4_ci_boot_perc <- ci_calc_boot_perc(d_sb_4_ci_boot, char_scale_d)
d_sb_8_ci_boot_perc <- ci_calc_boot_perc(d_sb_8_ci_boot, char_scale_d)

a_wv_cs_ci_perc <- ci_calc_boot_perc(a_wv_cs_ci, char_scale_a)
b_wv_cs_ci_perc <- ci_calc_boot_perc(b_wv_cs_ci, char_scale_b)
c_wv_cs_ci_perc <- ci_calc_boot_perc(c_wv_cs_ci, char_scale_c)
d_wv_cs_ci_perc <- ci_calc_boot_perc(d_wv_cs_ci, char_scale_d)

results <- rbind(
  a_wv_cs_ci_perc,
  a_sb_2_ci_boot_perc,
  a_sb_4_ci_boot_perc,
  a_sb_8_ci_boot_perc,
  b_wv_cs_ci_perc,
  b_sb_2_ci_boot_perc,
  b_sb_4_ci_boot_perc,
  b_sb_8_ci_boot_perc,
  c_wv_cs_ci_perc,
  c_sb_2_ci_boot_perc,
  c_sb_4_ci_boot_perc,
  c_sb_8_ci_boot_perc,
  d_wv_cs_ci_perc,
  d_sb_2_ci_boot_perc,
  d_sb_4_ci_boot_perc,
  d_sb_8_ci_boot_perc
)


for (i in 1:nrow(results)) {
  str_temp <- ""

  if ((i %% 4) == 1) {
    str_temp <- c("(a)", "(b)", "(c)", "(d)")[1 + (i %/% 4)]
  }

  cat(paste(
    c(str_temp, c(
      "Apr. Normal", "$(2\\log_2(n))^{-1}$ ", "$(4\\log_2(n))^{-1}$ ",
      "$(8\\log_2(n))^{-1}$ "
    )[((i - 1) %% 4) + 1], results[i, ]),
    collapse = " & "
  ))

  if (i < nrow(results)) {
    cat(" \\\\")
  }

  if ((i %% 4) == 0 && i != 16) {
    cat(" \\hline")
  }

  cat("\n")
}

################################################################################

# Auxiliary plot

set.seed(1)

ya <- model_a_sim(2048)
yb <- model_b_sim(2048)
yc <- model_c_sim(2048)
yd <- model_d_sim(2048)

wv_est_a <- wv_estimates(ya)
wv_est_b <- wv_estimates(yb)
wv_est_c <- wv_estimates(yc)
wv_est_d <- wv_estimates(yd)


{
  png(file = file.path(output_path, "Realizations.png"), width = 1200, height = 800, res = 150)

  par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 2.1))

  plot(log2(wv_est_a),
    xlab = "", ylab = "Var. de ondaletas (log)", pch = 19,
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xaxt = "n"
  )
  axis(1, labels = FALSE)

  abline(
    v = char_scale_a[[3]], lty = 2,
    lwd = 1.5
  )

  x <- 5:7
  y <- log2(wv_est_a)[5:7]
  quadratic <- solve(cbind(1, x, x^2), y)

  quadratic_fun <- function(x) {
    quadratic[[1]] + quadratic[[2]] * x + quadratic[[3]] * x^2
  }

  lines(seq(5, 7, length.out = 50), quadratic_fun(seq(5, 7, length.out = 50)),
    lwd = 1.5
  )

  plot(log2(wv_est_b),
    xlab = "", ylab = "", pch = 19,
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xaxt = "n"
  )
  axis(1, labels = FALSE)

  abline(
    v = char_scale_b[[3]], lty = 2,
    lwd = 1.5
  )

  x <- 3:5
  y <- log2(wv_est_b)[3:5]
  quadratic <- solve(cbind(1, x, x^2), y)

  lines(seq(3, 5, length.out = 50), quadratic_fun(seq(3, 5, length.out = 50)),
    lwd = 1.5
  )

  plot(log2(wv_est_c),
    xlab = "Nível", ylab = "Var. de ondaletas (log)", pch = 19,
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )

  abline(
    v = char_scale_c[[3]], lty = 2,
    lwd = 1.5
  )

  x <- 4:6
  y <- log2(wv_est_c[4:6])
  quadratic <- solve(cbind(1, x, x^2), y)

  lines(seq(4, 6, length.out = 50), quadratic_fun(seq(4, 6, length.out = 50)),
    lwd = 1.5
  )

  plot(log2(wv_est_d),
    xlab = "Nível", ylab = "", pch = 19,
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )

  dev.off()
}

save.image(file.path(workspace_dir, "6_char_scales.RData"))

################################################################################
