# =============================================================================
# 5_CI_wv.R
# =============================================================================
# Purpose  : Simulation study of bootstrap CI coverage for wavelet variance.
# Chapter  : Chapter 3, Appendix b.2
# Inputs   : None
# Outputs  : LaTeX tables for CI coverage.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "1_Simulation_functions.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------
source(file.path(base_path, "src", "2_Bootstrap_methods.R")) ################################################################################

# study for CI

set.seed(44)

b <- if (test_mode) 5 else 100
iterations <- if (test_mode) 2 else 100

a_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
b_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
c_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
d_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

a_wv_waveslim_gaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
b_wv_waveslim_gaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
c_wv_waveslim_gaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
d_wv_waveslim_gaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_waveslim_nongaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
b_wv_waveslim_nongaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
c_wv_waveslim_nongaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
d_wv_waveslim_nongaussian <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_waveslim_eta3 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
b_wv_waveslim_eta3 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
c_wv_waveslim_eta3 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)
d_wv_waveslim_eta3 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

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


for (iter in 1:iterations) {
  message(iter / iterations, "\r", appendLF = FALSE)
  flush.console()
  Sys.sleep(0.1)

  for (i in 1:3) {
    # print(c(iter,i))

    n <- 2^((2 * i - 1) + 6)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    # simulating from each model
    ya <- model_a_sim(n)
    yb <- model_b_sim(n)
    yc <- model_c_sim(n)
    yd <- model_d_sim(n)

    # calculating the point estimates of the wavelet variances
    a_wv[[i]] <- rbind(a_wv[[i]], wv_estimates(ya))
    b_wv[[i]] <- rbind(b_wv[[i]], wv_estimates(yb))
    c_wv[[i]] <- rbind(c_wv[[i]], wv_estimates(yc))
    d_wv[[i]] <- rbind(d_wv[[i]], wv_estimates(yd))

    # calculating the confidence intervals via waveslim
    a_wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(ya, n.levels = n_levels),
      p = 0.025, type = "gaussian"
    )[1:n_levels, ]
    b_wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(yb, n.levels = n_levels),
      p = 0.025, type = "gaussian"
    )[1:n_levels, ]
    c_wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(yc, n.levels = n_levels),
      p = 0.025, type = "gaussian"
    )[1:n_levels, ]
    d_wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(yd, n.levels = n_levels),
      p = 0.025, type = "gaussian"
    )[1:n_levels, ]

    a_wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(ya, n.levels = n_levels),
      p = 0.025, type = "nongaussian"
    )[1:n_levels, ]
    b_wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(yb, n.levels = n_levels),
      p = 0.025, type = "nongaussian"
    )[1:n_levels, ]
    c_wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(yc, n.levels = n_levels),
      p = 0.025, type = "nongaussian"
    )[1:n_levels, ]
    d_wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(yd, n.levels = n_levels),
      p = 0.025, type = "nongaussian"
    )[1:n_levels, ]


    a_wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(ya, n.levels = n_levels),
      p = 0.025, type = "eta3"
    )[1:n_levels, ]
    b_wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(yb, n.levels = n_levels),
      p = 0.025, type = "eta3"
    )[1:n_levels, ]
    c_wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(yc, n.levels = n_levels),
      p = 0.025, type = "eta3"
    )[1:n_levels, ]
    d_wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(yd, n.levels = n_levels),
      p = 0.025, type = "eta3"
    )[1:n_levels, ]

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

# confidence intervals

# ci_calculation_boot <- function(list1, alpha){
#
#   lower <- list("128" = NULL, "512" = NULL, "2048" = NULL)
#   upper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
#
#   for(iter in 1:iterations){
#
#     for(i in 1:3){
#
#       n <- 2^((2*i-1)+6)
#
#       n_levels <- floor(log2(1+(n-1)/(8-1)))
#
#       lower[[i]] <- rbind(lower[[i]], apply(list1[[i]][[iter]], 2, function(x){quantile(x, probs = alpha/2)}))
#
#       upper[[i]] <- rbind(upper[[i]], apply(list1[[i]][[iter]], 2, function(x){quantile(x, probs = 1 - (alpha/2))}))
#
#     }
#
#   }
#
#   return(list(lower, upper))
#
# }

#' Calculate bootstrap confidence intervals
#'
#' @param list1 List of bootstrap variances estimates
#' @param list2 List of point estimates for wavelet variances
#' @param alpha Significance level
#' @return A list containing lower and upper limits matrices
ci_calculation_boot <- function(list1, list2, alpha) {
  lower <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  upper <- list("128" = NULL, "512" = NULL, "2048" = NULL)

  for (iter in 1:iterations) {
    for (i in 1:3) {
      n <- 2^((2 * i - 1) + 6)

      n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

      # bootstrap variances estimates and point estimates for wavelet variances

      vars_boot <- apply(list1[[i]][[iter]], 2, function(x) {
        var(x)
      })

      mean_estimates <- list2[[i]][iter, ]

      # assuming that the wavelet variance estimators are proportional to a chi-squared
      # distribution, we may estimate these two parameters via

      a <- vars_boot / (2 * mean_estimates)

      eta <- (2 * mean_estimates^2) / vars_boot

      lower[[i]] <- rbind(lower[[i]], a * qchisq(rep(alpha / 2, n_levels), df = eta))

      upper[[i]] <- rbind(upper[[i]], a * qchisq(rep(1 - (alpha / 2), n_levels), df = eta))
    }
  }

  return(list(lower, upper))
}

a_wv_sb_2_ci <- ci_calculation_boot(a_wv_sb_2, a_wv, 0.05)
a_wv_sb_4_ci <- ci_calculation_boot(a_wv_sb_4, a_wv, 0.05)
a_wv_sb_8_ci <- ci_calculation_boot(a_wv_sb_8, a_wv, 0.05)

b_wv_sb_2_ci <- ci_calculation_boot(b_wv_sb_2, b_wv, 0.05)
b_wv_sb_4_ci <- ci_calculation_boot(b_wv_sb_4, b_wv, 0.05)
b_wv_sb_8_ci <- ci_calculation_boot(b_wv_sb_8, b_wv, 0.05)

c_wv_sb_2_ci <- ci_calculation_boot(c_wv_sb_2, c_wv, 0.05)
c_wv_sb_4_ci <- ci_calculation_boot(c_wv_sb_4, c_wv, 0.05)
c_wv_sb_8_ci <- ci_calculation_boot(c_wv_sb_8, c_wv, 0.05)

d_wv_sb_2_ci <- ci_calculation_boot(d_wv_sb_2, d_wv, 0.05)
d_wv_sb_4_ci <- ci_calculation_boot(d_wv_sb_4, d_wv, 0.05)
d_wv_sb_8_ci <- ci_calculation_boot(d_wv_sb_8, d_wv, 0.05)


#' Calculate true wavelet variances
#'
#' @param list1 Nested list of point estimates for wavelet variances
#' @return A vector of true wavelet variances
true_wv <- function(list1) {
  temp <- apply(list1[[3]], 2, function(x) {
    mean(x)
  })

  return(temp)
}

a_wv_mean <- true_wv(a_wv)
b_wv_mean <- true_wv(b_wv)
c_wv_mean <- true_wv(c_wv)
d_wv_mean <- true_wv(d_wv)

################################################################################

# coverages

#' Check CI coverage for bootstrap methods
#'
#' @param list1 List of logical coverage matrices
#' @param array2 Array of true values
#' @return List of coverage true/false
coverages_boot <- function(list1, array2) {
  output <- list(
    matrix(-1, nrow = iterations, ncol = 4),
    matrix(-1, nrow = iterations, ncol = 6),
    matrix(-1, nrow = iterations, ncol = 8)
  )

  for (i in 1:3) {
    n <- 2^((2 * i - 1) + 6)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    for (iter in 1:iterations) {
      for (j in 1:n_levels) {
        if (array2[j] >= list1[[1]][[i]][iter, j] & array2[j] <= list1[[2]][[i]][iter, j]) {
          output[[i]][iter, j] <- TRUE
        } else if (array2[j] < list1[[1]][[i]][iter, j] | array2[j] > list1[[2]][[i]][iter, j]) {
          output[[i]][iter, j] <- FALSE
        }
      }
    }
  }

  return(output)
}

a_wv_sb_2_ci_coverage <- coverages_boot(a_wv_sb_2_ci, a_wv_mean)
a_wv_sb_4_ci_coverage <- coverages_boot(a_wv_sb_4_ci, a_wv_mean)
a_wv_sb_8_ci_coverage <- coverages_boot(a_wv_sb_8_ci, a_wv_mean)

b_wv_sb_2_ci_coverage <- coverages_boot(b_wv_sb_2_ci, b_wv_mean)
b_wv_sb_4_ci_coverage <- coverages_boot(b_wv_sb_4_ci, b_wv_mean)
b_wv_sb_8_ci_coverage <- coverages_boot(b_wv_sb_8_ci, b_wv_mean)

c_wv_sb_2_ci_coverage <- coverages_boot(c_wv_sb_2_ci, c_wv_mean)
c_wv_sb_4_ci_coverage <- coverages_boot(c_wv_sb_4_ci, c_wv_mean)
c_wv_sb_8_ci_coverage <- coverages_boot(c_wv_sb_8_ci, c_wv_mean)

d_wv_sb_2_ci_coverage <- coverages_boot(d_wv_sb_2_ci, d_wv_mean)
d_wv_sb_4_ci_coverage <- coverages_boot(d_wv_sb_4_ci, d_wv_mean)
d_wv_sb_8_ci_coverage <- coverages_boot(d_wv_sb_8_ci, d_wv_mean)

#' Check CI coverage using waveslim
#'
#' @param list1 List of waveslim CI matrices
#' @param array2 Array of true values
#' @return List of coverage true/false
coverages_waveslim <- function(list1, array2) {
  output <- list(
    matrix(-1, nrow = iterations, ncol = 4),
    matrix(-1, nrow = iterations, ncol = 6),
    matrix(-1, nrow = iterations, ncol = 8)
  )

  for (i in 1:3) {
    n <- 2^((2 * i - 1) + 6)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    for (iter in 1:iterations) {
      for (j in 1:n_levels) {
        if (array2[j] >= list1[[i]][[iter]][j, 2] & array2[j] <= list1[[i]][[iter]][j, 3]) {
          output[[i]][iter, j] <- TRUE
        } else if (array2[j] < list1[[i]][[iter]][j, 2] | array2[j] > list1[[i]][[iter]][j, 2]) {
          output[[i]][iter, j] <- FALSE
        }
      }
    }
  }

  return(output)
}

a_wv_waveslim_gaussian_coverage <- coverages_waveslim(a_wv_waveslim_gaussian, a_wv_mean)
b_wv_waveslim_gaussian_coverage <- coverages_waveslim(b_wv_waveslim_gaussian, b_wv_mean)
c_wv_waveslim_gaussian_coverage <- coverages_waveslim(c_wv_waveslim_gaussian, c_wv_mean)
d_wv_waveslim_gaussian_coverage <- coverages_waveslim(d_wv_waveslim_gaussian, d_wv_mean)

a_wv_waveslim_nongaussian_coverage <- coverages_waveslim(a_wv_waveslim_nongaussian, a_wv_mean)
b_wv_waveslim_nongaussian_coverage <- coverages_waveslim(b_wv_waveslim_nongaussian, b_wv_mean)
c_wv_waveslim_nongaussian_coverage <- coverages_waveslim(c_wv_waveslim_nongaussian, c_wv_mean)
d_wv_waveslim_nongaussian_coverage <- coverages_waveslim(d_wv_waveslim_nongaussian, d_wv_mean)

a_wv_waveslim_eta3_coverage <- coverages_waveslim(a_wv_waveslim_eta3, a_wv_mean)
b_wv_waveslim_eta3_coverage <- coverages_waveslim(b_wv_waveslim_eta3, b_wv_mean)
c_wv_waveslim_eta3_coverage <- coverages_waveslim(c_wv_waveslim_eta3, c_wv_mean)
d_wv_waveslim_eta3_coverage <- coverages_waveslim(d_wv_waveslim_eta3, d_wv_mean)

################################################################################

# coverages (%)

#' Calculate coverage percentage
#'
#' @param list1 List of boolean coverages matrices
#' @return The coverage percentage vector
coverages_percentual <- function(list1) {
  temp <- list(NULL, NULL, NULL)

  for (i in 1:3) {
    temp[[i]] <- apply(list1[[i]], 2, function(x) {
      sum(x) / length(x)
    })
  }

  return(temp)
}

a_wv_sb_2_ci_coverage_percentual <- coverages_percentual(a_wv_sb_2_ci_coverage)
a_wv_sb_4_ci_coverage_percentual <- coverages_percentual(a_wv_sb_4_ci_coverage)
a_wv_sb_8_ci_coverage_percentual <- coverages_percentual(a_wv_sb_8_ci_coverage)

b_wv_sb_2_ci_coverage_percentual <- coverages_percentual(b_wv_sb_2_ci_coverage)
b_wv_sb_4_ci_coverage_percentual <- coverages_percentual(b_wv_sb_4_ci_coverage)
b_wv_sb_8_ci_coverage_percentual <- coverages_percentual(b_wv_sb_8_ci_coverage)

c_wv_sb_2_ci_coverage_percentual <- coverages_percentual(c_wv_sb_2_ci_coverage)
c_wv_sb_4_ci_coverage_percentual <- coverages_percentual(c_wv_sb_4_ci_coverage)
c_wv_sb_8_ci_coverage_percentual <- coverages_percentual(c_wv_sb_8_ci_coverage)

d_wv_sb_2_ci_coverage_percentual <- coverages_percentual(d_wv_sb_2_ci_coverage)
d_wv_sb_4_ci_coverage_percentual <- coverages_percentual(d_wv_sb_4_ci_coverage)
d_wv_sb_8_ci_coverage_percentual <- coverages_percentual(d_wv_sb_8_ci_coverage)

a_wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(a_wv_waveslim_gaussian_coverage)
b_wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(b_wv_waveslim_gaussian_coverage)
c_wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(c_wv_waveslim_gaussian_coverage)
d_wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(d_wv_waveslim_gaussian_coverage)

a_wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(a_wv_waveslim_nongaussian_coverage)
b_wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(b_wv_waveslim_nongaussian_coverage)
c_wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(c_wv_waveslim_nongaussian_coverage)
d_wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(d_wv_waveslim_nongaussian_coverage)

a_wv_waveslim_eta3_coverage_percentual <- coverages_percentual(a_wv_waveslim_eta3_coverage)
b_wv_waveslim_eta3_coverage_percentual <- coverages_percentual(b_wv_waveslim_eta3_coverage)
c_wv_waveslim_eta3_coverage_percentual <- coverages_percentual(c_wv_waveslim_eta3_coverage)
d_wv_waveslim_eta3_coverage_percentual <- coverages_percentual(d_wv_waveslim_eta3_coverage)

## ---- Tables and LateX Generation ----

#' Format simulation results to LaTeX rules
#'
#' @param results Aggregated matrix of percentages
#' @return Latex table row printed to console
latex_fun <- function(results) {
  for (i in 1:nrow(results)) {
    str_temp <- ""

    if ((i + 1) %% 3 == 0) {
      str_temp <- c(
        "Gaussian ", "$\\hat{\\eta}_3$ ", "Multitaper ", "$(2\\log_2(n))^{-1}$ ",
        "$(4\\log_2(n))^{-1}$ ", "$(8\\log_2(n))^{-1}$ "
      )[(i + 1) %/% 3]
    }

    cat(paste(c(str_temp, c(128, 512, 2048)[((i - 1) %% 3) + 1], results[i, ]), collapse = " & "))

    if (i < nrow(results)) {
      cat(" \\\\")
    }


    if ((i %% 3) == 0 && i != 18) {
      cat("\\hline")
    }

    cat("\n")
  }
}

results_a <- rbind(
  c(a_wv_waveslim_gaussian_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_waveslim_gaussian_coverage_percentual[[2]], rep("", 2)),
  a_wv_waveslim_gaussian_coverage_percentual[[3]],
  c(a_wv_waveslim_eta3_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_waveslim_eta3_coverage_percentual[[2]], rep("", 2)),
  a_wv_waveslim_eta3_coverage_percentual[[3]],
  c(a_wv_waveslim_nongaussian_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_waveslim_nongaussian_coverage_percentual[[2]], rep("", 2)),
  a_wv_waveslim_nongaussian_coverage_percentual[[3]],
  c(a_wv_sb_2_ci_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_sb_2_ci_coverage_percentual[[2]], rep("", 2)),
  a_wv_sb_2_ci_coverage_percentual[[3]],
  c(a_wv_sb_4_ci_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_sb_4_ci_coverage_percentual[[2]], rep("", 2)),
  a_wv_sb_4_ci_coverage_percentual[[3]],
  c(a_wv_sb_8_ci_coverage_percentual[[1]], rep("", 4)),
  c(a_wv_sb_8_ci_coverage_percentual[[2]], rep("", 2)),
  a_wv_sb_8_ci_coverage_percentual[[3]]
)

results_b <- rbind(
  c(b_wv_waveslim_gaussian_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_waveslim_gaussian_coverage_percentual[[2]], rep("", 2)),
  b_wv_waveslim_gaussian_coverage_percentual[[3]],
  c(b_wv_waveslim_eta3_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_waveslim_eta3_coverage_percentual[[2]], rep("", 2)),
  b_wv_waveslim_eta3_coverage_percentual[[3]],
  c(b_wv_waveslim_nongaussian_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_waveslim_nongaussian_coverage_percentual[[2]], rep("", 2)),
  b_wv_waveslim_nongaussian_coverage_percentual[[3]],
  c(b_wv_sb_2_ci_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_sb_2_ci_coverage_percentual[[2]], rep("", 2)),
  b_wv_sb_2_ci_coverage_percentual[[3]],
  c(b_wv_sb_4_ci_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_sb_4_ci_coverage_percentual[[2]], rep("", 2)),
  b_wv_sb_4_ci_coverage_percentual[[3]],
  c(b_wv_sb_8_ci_coverage_percentual[[1]], rep("", 4)),
  c(b_wv_sb_8_ci_coverage_percentual[[2]], rep("", 2)),
  b_wv_sb_8_ci_coverage_percentual[[3]]
)

results_c <- rbind(
  c(c_wv_waveslim_gaussian_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_waveslim_gaussian_coverage_percentual[[2]], rep("", 2)),
  c_wv_waveslim_gaussian_coverage_percentual[[3]],
  c(c_wv_waveslim_eta3_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_waveslim_eta3_coverage_percentual[[2]], rep("", 2)),
  c_wv_waveslim_eta3_coverage_percentual[[3]],
  c(c_wv_waveslim_nongaussian_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_waveslim_nongaussian_coverage_percentual[[2]], rep("", 2)),
  c_wv_waveslim_nongaussian_coverage_percentual[[3]],
  c(c_wv_sb_2_ci_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_sb_2_ci_coverage_percentual[[2]], rep("", 2)),
  c_wv_sb_2_ci_coverage_percentual[[3]],
  c(c_wv_sb_4_ci_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_sb_4_ci_coverage_percentual[[2]], rep("", 2)),
  c_wv_sb_4_ci_coverage_percentual[[3]],
  c(c_wv_sb_8_ci_coverage_percentual[[1]], rep("", 4)),
  c(c_wv_sb_8_ci_coverage_percentual[[2]], rep("", 2)),
  c_wv_sb_8_ci_coverage_percentual[[3]]
)

results_d <- rbind(
  c(d_wv_waveslim_gaussian_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_waveslim_gaussian_coverage_percentual[[2]], rep("", 2)),
  d_wv_waveslim_gaussian_coverage_percentual[[3]],
  c(d_wv_waveslim_eta3_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_waveslim_eta3_coverage_percentual[[2]], rep("", 2)),
  d_wv_waveslim_eta3_coverage_percentual[[3]],
  c(d_wv_waveslim_nongaussian_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_waveslim_nongaussian_coverage_percentual[[2]], rep("", 2)),
  d_wv_waveslim_nongaussian_coverage_percentual[[3]],
  c(d_wv_sb_2_ci_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_sb_2_ci_coverage_percentual[[2]], rep("", 2)),
  d_wv_sb_2_ci_coverage_percentual[[3]],
  c(d_wv_sb_4_ci_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_sb_4_ci_coverage_percentual[[2]], rep("", 2)),
  d_wv_sb_4_ci_coverage_percentual[[3]],
  c(d_wv_sb_8_ci_coverage_percentual[[1]], rep("", 4)),
  c(d_wv_sb_8_ci_coverage_percentual[[2]], rep("", 2)),
  d_wv_sb_8_ci_coverage_percentual[[3]]
)

latex_fun(results_a)
latex_fun(results_b)
latex_fun(results_c)
latex_fun(results_d)

################################################################################


#' Detailed LaTeX formatting wrapper function
#'
#' @param results_a Results from Model A
#' @param results_b Results from Model b
#' @param results_c Results from Model C
#' @param results_d Results from Model D
#' @return Latex table row median values printed to console
latex_fun2 <- function(results_a, results_b, results_c, results_d) {
  for (i in 1:nrow(results_a)) {
    if (i %% 3 == 1) {
      val1a <- median(as.numeric(results_a[i, 1:2]), na.rm = TRUE)
      val2a <- median(as.numeric(results_a[i, 3:4]), na.rm = TRUE)

      val1b <- median(as.numeric(results_b[i, 1:2]), na.rm = TRUE)
      val2b <- median(as.numeric(results_b[i, 3:4]), na.rm = TRUE)

      val1c <- median(as.numeric(results_c[i, 1:2]), na.rm = TRUE)
      val2c <- median(as.numeric(results_c[i, 3:4]), na.rm = TRUE)

      val1d <- median(as.numeric(results_d[i, 1:2]), na.rm = TRUE)
      val2d <- median(as.numeric(results_d[i, 3:4]), na.rm = TRUE)
    } else if (i %% 3 == 2) {
      val1a <- median(as.numeric(results_a[i, 1:3]), na.rm = TRUE)
      val2a <- median(as.numeric(results_a[i, 4:6]), na.rm = TRUE)

      val1b <- median(as.numeric(results_b[i, 1:3]), na.rm = TRUE)
      val2b <- median(as.numeric(results_b[i, 4:6]), na.rm = TRUE)

      val1c <- median(as.numeric(results_c[i, 1:3]), na.rm = TRUE)
      val2c <- median(as.numeric(results_c[i, 4:6]), na.rm = TRUE)

      val1d <- median(as.numeric(results_d[i, 1:3]), na.rm = TRUE)
      val2d <- median(as.numeric(results_d[i, 4:6]), na.rm = TRUE)
    } else if (i %% 3 == 0) {
      val1a <- median(as.numeric(results_a[i, 1:4]), na.rm = TRUE)
      val2a <- median(as.numeric(results_a[i, 5:8]), na.rm = TRUE)

      val1b <- median(as.numeric(results_b[i, 1:4]), na.rm = TRUE)
      val2b <- median(as.numeric(results_b[i, 5:8]), na.rm = TRUE)

      val1c <- median(as.numeric(results_c[i, 1:4]), na.rm = TRUE)
      val2c <- median(as.numeric(results_c[i, 5:8]), na.rm = TRUE)

      val1d <- median(as.numeric(results_d[i, 1:4]), na.rm = TRUE)
      val2d <- median(as.numeric(results_d[i, 5:8]), na.rm = TRUE)
    }

    vals <- c(val1a, val2a, val1b, val2b, val1c, val2c, val1d, val2d)


    str_temp <- ""

    if ((i + 1) %% 3 == 0) {
      str_temp <- c(
        "Gaussian ", "$\\hat{\\eta}_3$ ", "Multitaper ", "$(2\\log_2(n))^{-1}$ ",
        "$(4\\log_2(n))^{-1}$ ", "$(8\\log_2(n))^{-1}$ "
      )[(i + 1) %/% 3]
    }

    cat(paste(c(str_temp, c(128, 512, 2048)[((i - 1) %% 3) + 1], vals), collapse = " & "))

    if (i < nrow(results_a)) {
      cat(" \\\\")
    }


    if ((i %% 3) == 0 && i != 18) {
      cat("\\hline")
    }

    cat("\n")
  }
}

latex_fun2(results_a, results_b, results_c, results_d)

save.image(file.path(workspace_dir, "5_CI_wv.RData"))
