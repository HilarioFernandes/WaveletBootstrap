# =============================================================================
# 9_Comparison_sim.R
# =============================================================================
# Purpose  : Simulation study for wavelet variance ratio tests (F-test vs Bootstrap).
# Chapter  : Chapter 3, Appendix b.4
# Inputs   : None
# Outputs  : Size and power rejection proportions for ratio tests.
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

alpha_0 <- 0.05

#' Run simulation study comparing F-test and Bootstrap for wavelet variance ratios
#'
#' @param n Sample size
#' @param mult_fact Variance multiplier for Group y
#' @param indep_ground Logical. If TRUE, groups are independent processes.
#' @param indep_method Logical. If TRUE, bootstrap resamples groups independently.
#' @param model Model identifier ("A" through "D")
#' @param iterations Number of simulation iterations
#' @param b Number of bootstrap reps
#' @param alpha_0 Target significance level
#' @return A list containing rejection proportions for F-test and Bootstrap
simulations_comparison <- function(n, mult_fact, indep_ground, indep_method, model, iterations, b, alpha_0) {
  rejections_boot <- NULL

  rejections_f <- NULL

  for (iter in 1:iterations) {
    print(iter)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    sim_fun <- get(paste0("model_", model, "_sim"))

    if (indep_ground == TRUE) {
      x <- sim_fun(n)
      y <- mult_fact * sim_fun(n)
    } else {
      x <- sim_fun(2 * n)
      y <- mult_fact * x[(n + 1):(2 * n)]
      x <- x[1:n]
    }

    x_wv <- wv_estimates(x)
    y_wv <- wv_estimates(y)

    ratios <- y_wv / x_wv

    x_y <- c(x, y)

    # F approximation-based tests

    ratios_f <- NULL

    rejections_f_temp <- NULL

    for (jmax in 1:n_levels) {
      alpha <- alpha_0 / jmax

      js <- 1:jmax

      l_js <- (2^(js) - 1) * (8 - 1) + 1

      n_js <- n - l_js + 1

      eta_js <- sapply(n_js / 2^js, function(x) {
        max(x, 1)
      })

      quants_f <- t(sapply(eta_js, function(x) {
        qf(c(alpha / 2, 1 - alpha / 2), x, x)
      }))

      rejections_f_temp <- c(rejections_f_temp, max(apply(cbind(ratios[1:jmax], quants_f), 1, interval_check)))
    }

    rejections_f <- rbind(rejections_f, rejections_f_temp)

    rejection_f_proportions <- apply(rejections_f, 2, function(x) {
      sum(x) / length(x)
    })

    # bootstrap-based tests

    ratios_boot <- NULL

    for (b in 1:b) {
      if (indep_method == FALSE) {
        boot_indexes <- block_boot(2 * n, "SB", 1 / (4 * log2(2 * n)))

        x_y_boot <- x_y[boot_indexes]
        x_boot_wv <- wv_estimates(x_y_boot[1:n])
        y_boot_wv <- wv_estimates(x_y_boot[(n + 1):(2 * n)])
      } else {
        boot_indexes_a <- block_boot(n, "SB", 1 / (4 * log2(n)))
        boot_indexes_b <- block_boot(n, "SB", 1 / (4 * log2(n)))

        x_boot <- x[boot_indexes_a]
        y_boot <- y[boot_indexes_b]

        x_boot_wv <- wv_estimates(x_boot)
        y_boot_wv <- wv_estimates(y_boot)
      }

      ratios_boot <- rbind(ratios_boot, array(y_boot_wv / x_boot_wv))
    }

    rejections_boot_temp <- NULL

    for (jmax in 1:n_levels) {
      alpha <- alpha_0 / jmax

      quants_boot <- t(apply(ratios_boot, 2, function(x) {
        quantile(x, c(alpha / 2, 1 - alpha / 2))
      })[, 1:jmax])

      rejections_boot_temp <- c(rejections_boot_temp, max(apply(cbind(ratios[1:jmax], quants_boot), 1, interval_check)))
    }

    rejections_boot <- rbind(rejections_boot, rejections_boot_temp)

    rejection_boot_proportions <- apply(rejections_boot, 2, function(x) {
      sum(x) / length(x)
    })
  }

  return(list(rejection_f_proportions, rejection_boot_proportions))
}

# quick tests
if (FALSE) {
  # simulations_comparison(128, 1, TRUE, TRUE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, TRUE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, TRUE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, TRUE, "D", 20,20, alpha_0)
  #
  # #mau comportamento
  # simulations_comparison(128, 1.5, TRUE, TRUE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, TRUE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, TRUE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, TRUE, "D", 20,20, alpha_0)
  #
  # #bom comportamento
  # simulations_comparison(128, 1, TRUE, FALSE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, FALSE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, FALSE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1, TRUE, FALSE, "D", 20,20, alpha_0)
  #
  # #bom comportamento
  # simulations_comparison(128, 1.5, TRUE, FALSE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, FALSE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, FALSE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, TRUE, FALSE, "D", 20,20, alpha_0)
  #
  # #bom comportamento
  # simulations_comparison(128, 1, FALSE, FALSE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, FALSE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, FALSE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, FALSE, "D", 20,20, alpha_0)
  #
  # #bom comportamento
  # simulations_comparison(128, 1.5, FALSE, FALSE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, FALSE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, FALSE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, FALSE, "D", 20,20, alpha_0)
  #
  # #bom comportamento
  # simulations_comparison(128, 1, FALSE, TRUE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, TRUE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, TRUE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1, FALSE, TRUE, "D", 20,20, alpha_0)
  #
  # #mau comportamento
  # simulations_comparison(128, 1.5, FALSE, TRUE, "A", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, TRUE, "B", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, TRUE, "C", 20,20, alpha_0)
  # simulations_comparison(128, 1.5, FALSE, TRUE, "D", 20,20, alpha_0)
}

################################################################################

iterations <- if (test_mode) 2 else 100
b <- if (test_mode) 2 else 100

################################################################################

set.seed(58)

# No variance change, independent blocks, independent method

start_time <- Sys.time()

# model A

rejection_proportions_f_a_1_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_A_1_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, TRUE, "A", iterations, b, alpha_0)

  rejection_proportions_f_a_1_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_A_1_indepground_indepmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_f_b_1_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_B_1_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, TRUE, "B", iterations, b, alpha_0)

  rejection_proportions_f_b_1_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_B_1_indepground_indepmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_f_c_1_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_C_1_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, TRUE, "C", iterations, b, alpha_0)

  rejection_proportions_f_c_1_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_C_1_indepground_indepmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_f_d_1_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_D_1_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, TRUE, "D", iterations, b, alpha_0)

  rejection_proportions_f_d_1_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_D_1_indepground_indepmethod[[i]] <- temp[[2]]
}

  end_time <- Sys.time()

  duration <- end_time - start_time
  print(duration)

save.image(file.path(workspace_dir, "9_Comparison_sim_part_1.RData"))
set.seed(59)

# variance change (factor 1.5), independent blocks, independent method

start_time <- Sys.time()

# model A

rejection_proportions_f_a_1_5_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_a_1_5_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, TRUE, "A", iterations, b, alpha_0)

  rejection_proportions_f_a_1_5_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_a_1_5_indepground_indepmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_f_b_1_5_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_b_1_5_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, TRUE, "B", iterations, b, alpha_0)

  rejection_proportions_f_b_1_5_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_b_1_5_indepground_indepmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_f_c_1_5_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_c_1_5_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, TRUE, "C", iterations, b, alpha_0)

  rejection_proportions_f_c_1_5_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_c_1_5_indepground_indepmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_f_d_1_5_indepground <- list(NULL, NULL, NULL)
rejection_proportions_boot_d_1_5_indepground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, TRUE, "D", iterations, b, alpha_0)

  rejection_proportions_f_d_1_5_indepground[[i]] <- temp[[1]]
  rejection_proportions_boot_d_1_5_indepground_indepmethod[[i]] <- temp[[2]]
}

  end_time <- Sys.time()

  duration <- end_time - start_time
  print(duration)

################################################################################

save.image(file.path(workspace_dir, "9_Comparison_sim_part_2.RData"))
set.seed(60)

# No variance change, independent blocks, dependent method

start_time <- Sys.time()

# model A

rejection_proportions_boot_A_1_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, FALSE, "A", iterations, b, alpha_0)

  rejection_proportions_boot_A_1_indepground_depmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_boot_B_1_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, FALSE, "B", iterations, b, alpha_0)

  rejection_proportions_boot_B_1_indepground_depmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_boot_C_1_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, FALSE, "C", iterations, b, alpha_0)

  rejection_proportions_boot_C_1_indepground_depmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_boot_D_1_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, TRUE, FALSE, "D", iterations, b, alpha_0)

  rejection_proportions_boot_D_1_indepground_depmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "9_Comparison_sim_part_3.RData"))
set.seed(61)

# variance change (factor 1.5), independent blocks, dependent method

start_time <- Sys.time()

# model A

rejection_proportions_boot_a_1_5_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, FALSE, "A", iterations, b, alpha_0)

  rejection_proportions_boot_a_1_5_indepground_depmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_boot_b_1_5_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, FALSE, "B", iterations, b, alpha_0)

  rejection_proportions_boot_b_1_5_indepground_depmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_boot_c_1_5_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, FALSE, "C", iterations, b, alpha_0)

  rejection_proportions_boot_c_1_5_indepground_depmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_boot_d_1_5_indepground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, TRUE, FALSE, "D", iterations, b, alpha_0)

  rejection_proportions_boot_d_1_5_indepground_depmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

################################################################################

save.image(file.path(workspace_dir, "9_Comparison_sim_part_4.RData"))
set.seed(62)

# No variance change, dependent blocks, dependent method

start_time <- Sys.time()

# model A

rejection_proportions_f_a_1_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_A_1_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, FALSE, "A", iterations, b, alpha_0)

  rejection_proportions_f_a_1_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_A_1_depground_depmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_f_b_1_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_B_1_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, FALSE, "B", iterations, b, alpha_0)

  rejection_proportions_f_b_1_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_B_1_depground_depmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_f_c_1_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_C_1_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, FALSE, "C", iterations, b, alpha_0)

  rejection_proportions_f_c_1_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_C_1_depground_depmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_f_d_1_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_D_1_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, FALSE, "D", iterations, b, alpha_0)

  rejection_proportions_f_d_1_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_D_1_depground_depmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "9_Comparison_sim_part_5.RData"))
set.seed(63)

# variance change (factor 1.5), dependent blocks, dependent method

start_time <- Sys.time()

# model A

rejection_proportions_f_a_1_5_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_a_1_5_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, FALSE, "A", iterations, b, alpha_0)

  rejection_proportions_f_a_1_5_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_a_1_5_depground_depmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_f_b_1_5_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_b_1_5_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, FALSE, "B", iterations, b, alpha_0)

  rejection_proportions_f_b_1_5_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_b_1_5_depground_depmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_f_c_1_5_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_c_1_5_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, FALSE, "C", iterations, b, alpha_0)

  rejection_proportions_f_c_1_5_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_c_1_5_depground_depmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_f_d_1_5_depground <- list(NULL, NULL, NULL)
rejection_proportions_boot_d_1_5_depground_depmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, FALSE, "D", iterations, b, alpha_0)

  rejection_proportions_f_d_1_5_depground[[i]] <- temp[[1]]
  rejection_proportions_boot_d_1_5_depground_depmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

################################################################################

save.image(file.path(workspace_dir, "9_Comparison_sim_part_6.RData"))
set.seed(64)

# No variance change, dependent blocks, independent method

start_time <- Sys.time()

# model A

rejection_proportions_boot_A_1_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, TRUE, "A", iterations, b, alpha_0)

  rejection_proportions_boot_A_1_depground_indepmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_boot_B_1_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, TRUE, "B", iterations, b, alpha_0)

  rejection_proportions_boot_B_1_depground_indepmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_boot_C_1_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, TRUE, "C", iterations, b, alpha_0)

  rejection_proportions_boot_C_1_depground_indepmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_boot_D_1_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1, FALSE, TRUE, "D", iterations, b, alpha_0)

  rejection_proportions_boot_D_1_depground_indepmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "9_Comparison_sim_part_7.RData"))
set.seed(65)

# variance change (factor 1.5), dependent blocks, independent method

start_time <- Sys.time()

# model A

rejection_proportions_boot_a_1_5_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, TRUE, "A", iterations, b, alpha_0)

  rejection_proportions_boot_a_1_5_depground_indepmethod[[i]] <- temp[[2]]
}

# model b

rejection_proportions_boot_b_1_5_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, TRUE, "B", iterations, b, alpha_0)

  rejection_proportions_boot_b_1_5_depground_indepmethod[[i]] <- temp[[2]]
}

# model C

rejection_proportions_boot_c_1_5_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, TRUE, "C", iterations, b, alpha_0)

  rejection_proportions_boot_c_1_5_depground_indepmethod[[i]] <- temp[[2]]
}

# model D

rejection_proportions_boot_d_1_5_depground_indepmethod <- list(NULL, NULL, NULL)

for (i in 1:3) {
  n <- 2^((2 * i - 1) + 6)

  temp <- simulations_comparison(n, 1.5, FALSE, TRUE, "D", iterations, b, alpha_0)

  rejection_proportions_boot_d_1_5_depground_indepmethod[[i]] <- temp[[2]]
}

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

################################################################################

# bom comportamento
rejection_proportions_boot_A_1_indepground_indepmethod
rejection_proportions_boot_B_1_indepground_indepmethod
rejection_proportions_boot_C_1_indepground_indepmethod
rejection_proportions_boot_D_1_indepground_indepmethod

# mau comportamento
rejection_proportions_boot_a_1_5_indepground_indepmethod
rejection_proportions_boot_b_1_5_indepground_indepmethod
rejection_proportions_boot_c_1_5_indepground_indepmethod
rejection_proportions_boot_d_1_5_indepground_indepmethod

# bom comportamento
rejection_proportions_boot_A_1_indepground_depmethod
rejection_proportions_boot_B_1_indepground_depmethod
rejection_proportions_boot_C_1_indepground_depmethod
rejection_proportions_boot_D_1_indepground_depmethod

# bom comportamento
rejection_proportions_boot_a_1_5_indepground_depmethod
rejection_proportions_boot_b_1_5_indepground_depmethod
rejection_proportions_boot_c_1_5_indepground_depmethod
rejection_proportions_boot_d_1_5_indepground_depmethod

# bom comportamento
rejection_proportions_boot_A_1_depground_depmethod
rejection_proportions_boot_B_1_depground_depmethod
rejection_proportions_boot_C_1_depground_depmethod
rejection_proportions_boot_D_1_depground_depmethod

# bom comportamento
rejection_proportions_boot_a_1_5_depground_depmethod
rejection_proportions_boot_b_1_5_depground_depmethod
rejection_proportions_boot_c_1_5_depground_depmethod
rejection_proportions_boot_d_1_5_depground_depmethod

# bom comportamento
rejection_proportions_boot_A_1_depground_indepmethod
rejection_proportions_boot_B_1_depground_indepmethod
rejection_proportions_boot_C_1_depground_indepmethod
rejection_proportions_boot_D_1_depground_indepmethod

# mau comportamento
rejection_proportions_boot_a_1_5_depground_indepmethod
rejection_proportions_boot_b_1_5_depground_indepmethod
rejection_proportions_boot_c_1_5_depground_indepmethod
rejection_proportions_boot_d_1_5_depground_indepmethod

# F approximation
rejection_proportions_f_a_1_indepground
rejection_proportions_f_b_1_indepground
rejection_proportions_f_c_1_indepground
rejection_proportions_f_d_1_indepground

rejection_proportions_f_a_1_5_indepground
rejection_proportions_f_b_1_5_indepground
rejection_proportions_f_c_1_5_indepground
rejection_proportions_f_d_1_5_indepground

rejection_proportions_f_a_1_depground
rejection_proportions_f_b_1_depground
rejection_proportions_f_c_1_depground
rejection_proportions_f_d_1_depground

rejection_proportions_f_a_1_5_depground
rejection_proportions_f_b_1_5_depground
rejection_proportions_f_c_1_5_depground
rejection_proportions_f_d_1_5_depground

################################################################################

#' Generate detailed LaTeX table for rejection rates
#'
#' @param list1,list2,list3,list4 Rejection proportions (size) for models A, b, C, D
#' @param list5,list6,list7,list8 Rejection proportions (power) for models A, b, C, D
#' @return Prints a formatted LaTeX table to the console
latex_fun <- function(list1, list2, list3, list4, list5, list6, list7, list8) {
  cat("\\begin{table}[h!] \n")

  cat("\\centering \n")

  cat("\\caption{Taxas empíricas de rejeição de $H_0$ (tamanho/poder fora/dentro de parênteses) para blocos \\textcolor{red}{independentes/dependentes} e utilizando o \\textit{bootstrap} \\textcolor{red}{com/sem} discriminação de blocos.} \n")

  cat("\\label{tabela_compar} \n")

  cat("\\begin{tabular}{lcccccccc} \n")

  cat("\\multicolumn{1}{l|}{Modelo (a)} & $\\# j = 1$ & $\\# j = 2$ & $\\# j = 3$ & $\\# j = 4$ & $\\# j = 5$ & $\\# j = 6$ & $\\# j = 7$ & $\\# j = 8$ \\\\ \\hline \n")

  for (i in 1:length(list1)) {
    if (i == 1) {
      cat(paste0("\\multicolumn{1}{l|}{$n=128$} & ", paste(list1[[i]], collapse = " & "), " & & & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list5[[i]], collapse = ") & ("), ") & & & & \\\\ \n"))
    } else if (i == 2) {
      cat(paste0("\\multicolumn{1}{l|}{$n=512$} & ", paste(list1[[i]], collapse = " & "), " & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list5[[i]], collapse = ") & ("), ") & & \\\\ \n"))
    } else if (i == 3) {
      cat(paste0("\\multicolumn{1}{l|}{$n=2048$} & ", paste(list1[[i]], collapse = " & "), " \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list5[[i]], collapse = ") & ("), ") \\\\ \n"))
    }
  }

  cat("& & & & & & & & \\\\ \n")

  cat("\\multicolumn{1}{l|}{Modelo (b)} & & & & & & & & \\\\ \\hline \n")

  for (i in 1:length(list2)) {
    if (i == 1) {
      cat(paste0("\\multicolumn{1}{l|}{$n=128$} & ", paste(list2[[i]], collapse = " & "), " & & & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list6[[i]], collapse = ") & ("), ") & & & & \\\\ \n"))
    } else if (i == 2) {
      cat(paste0("\\multicolumn{1}{l|}{$n=512$} & ", paste(list2[[i]], collapse = " & "), " & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list2[[i]], collapse = ") & ("), ") & & \\\\ \n"))
    } else if (i == 3) {
      cat(paste0("\\multicolumn{1}{l|}{$n=2048$} & ", paste(list2[[i]], collapse = " & "), " \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list6[[i]], collapse = ") & ("), ") \\\\ \n"))
    }
  }

  cat("& & & & & & & & \\\\ \n")

  cat("\\multicolumn{1}{l|}{Modelo (c)} & & & & & & & & \\\\ \\hline \n")

  for (i in 1:length(list3)) {
    if (i == 1) {
      cat(paste0("\\multicolumn{1}{l|}{$n=128$} & ", paste(list3[[i]], collapse = " & "), " & & & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list7[[i]], collapse = ") & ("), ") & & & & \\\\ \n"))
    } else if (i == 2) {
      cat(paste0("\\multicolumn{1}{l|}{$n=512$} & ", paste(list3[[i]], collapse = " & "), " & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list7[[i]], collapse = ") & ("), ") & & \\\\ \n"))
    } else if (i == 3) {
      cat(paste0("\\multicolumn{1}{l|}{$n=2048$} & ", paste(list3[[i]], collapse = " & "), " \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list7[[i]], collapse = ") & ("), ") \\\\ \n"))
    }
  }

  cat("& & & & & & & & \\\\ \n")

  cat("\\multicolumn{1}{l|}{Modelo (d)} & & & & & & & & \\\\ \\hline \n")

  for (i in 1:length(list4)) {
    if (i == 1) {
      cat(paste0("\\multicolumn{1}{l|}{$n=128$} & ", paste(list4[[i]], collapse = " & "), " & & & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list8[[i]], collapse = ") & ("), ") & & & & \\\\ \n"))
    } else if (i == 2) {
      cat(paste0("\\multicolumn{1}{l|}{$n=512$} & ", paste(list4[[i]], collapse = " & "), " & & \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list8[[i]], collapse = ") & ("), ") & & \\\\ \n"))
    } else if (i == 3) {
      cat(paste0("\\multicolumn{1}{l|}{$n=2048$} & ", paste(list4[[i]], collapse = " & "), " \\\\ \n"))

      cat(paste0("\\multicolumn{1}{l|}{} & (", paste(list8[[i]], collapse = ") & ("), ") \n"))
    }
  }

  cat("\\end{tabular} \n")

  cat("\\end{table} \n")
}

latex_fun(
  rejection_proportions_boot_A_1_indepground_indepmethod,
  rejection_proportions_boot_B_1_indepground_indepmethod,
  rejection_proportions_boot_C_1_indepground_indepmethod,
  rejection_proportions_boot_D_1_indepground_indepmethod,
  rejection_proportions_boot_a_1_5_indepground_indepmethod,
  rejection_proportions_boot_b_1_5_indepground_indepmethod,
  rejection_proportions_boot_c_1_5_indepground_indepmethod,
  rejection_proportions_boot_d_1_5_indepground_indepmethod
)

latex_fun(
  rejection_proportions_boot_A_1_depground_indepmethod,
  rejection_proportions_boot_B_1_depground_indepmethod,
  rejection_proportions_boot_C_1_depground_indepmethod,
  rejection_proportions_boot_D_1_depground_indepmethod,
  rejection_proportions_boot_a_1_5_depground_indepmethod,
  rejection_proportions_boot_b_1_5_depground_indepmethod,
  rejection_proportions_boot_c_1_5_depground_indepmethod,
  rejection_proportions_boot_d_1_5_depground_indepmethod
)

latex_fun(
  rejection_proportions_boot_A_1_indepground_depmethod,
  rejection_proportions_boot_B_1_indepground_depmethod,
  rejection_proportions_boot_C_1_indepground_depmethod,
  rejection_proportions_boot_D_1_indepground_depmethod,
  rejection_proportions_boot_a_1_5_indepground_depmethod,
  rejection_proportions_boot_b_1_5_indepground_depmethod,
  rejection_proportions_boot_c_1_5_indepground_depmethod,
  rejection_proportions_boot_d_1_5_indepground_depmethod
)

latex_fun(
  rejection_proportions_boot_A_1_depground_depmethod,
  rejection_proportions_boot_B_1_depground_depmethod,
  rejection_proportions_boot_C_1_depground_depmethod,
  rejection_proportions_boot_D_1_depground_depmethod,
  rejection_proportions_boot_a_1_5_depground_depmethod,
  rejection_proportions_boot_b_1_5_depground_depmethod,
  rejection_proportions_boot_c_1_5_depground_depmethod,
  rejection_proportions_boot_d_1_5_depground_depmethod
)

latex_fun(
  rejection_proportions_f_a_1_indepground,
  rejection_proportions_f_b_1_indepground,
  rejection_proportions_f_c_1_indepground,
  rejection_proportions_f_d_1_indepground,
  rejection_proportions_f_a_1_5_indepground,
  rejection_proportions_f_b_1_5_indepground,
  rejection_proportions_f_c_1_5_indepground,
  rejection_proportions_f_d_1_5_indepground
)

latex_fun(
  rejection_proportions_f_a_1_depground,
  rejection_proportions_f_b_1_depground,
  rejection_proportions_f_c_1_depground,
  rejection_proportions_f_d_1_depground,
  rejection_proportions_f_a_1_5_depground,
  rejection_proportions_f_b_1_5_depground,
  rejection_proportions_f_c_1_5_depground,
  rejection_proportions_f_d_1_5_depground
)

################################################################################

# Table with rule of thumb for level choice

#' Generate summary LaTeX table for specific level choice (rule of thumb)
#'
#' @param list1,list2,list3,list4 Bootstrap rejection rates (H0) for 4 models
#' @param list5,list6,list7,list8 Bootstrap rejection rates (H1) for 4 models
#' @param list9,list10,list11,list12 F-approximation rates (H0) for 4 models
#' @param list13,list14,list15,list16 F-approximation rates (H1) for 4 models
#' @return Prints a summary LaTeX table comparing methods to the console
latex_fun2 <- function(list1, list2, list3, list4, list5, list6, list7, list8,
                       list9, list10, list11, list12, list13, list14, list15, list16) {
  cat("\\begin{table}[h!] \n")

  cat("\\centering \n")

  cat("\\caption{Taxas empíricas de rejeição de $H_0$ (poder indicado dentro de parênteses) para blocos \\textcolor{red}{independentes/dependentes}.} \n")

  cat("\\label{tabela_compar} \n")

  cat("\\begin{tabular}{l|llllllll} \n")

  cat("Método   & Boot. & Aprox.$F$ & Boot. & Aprox.$F$ & Boot. & Aprox.$F$ & Boot. & Aprox.$F$ \\\\ \n")

  cat("Modelo   & (a)   & (a)       & (b)   & (b)       & (c)   & (c)       & (d)   & (d)       \\\\ \\hline \n")

  cat(paste0(
    "$n=128$  & ", list1[[1]][2], " & ", list9[[1]][2], " & ",
    list2[[1]][2], " & ", list10[[1]][2], " & ",
    list3[[1]][2], " & ", list11[[1]][2], " & ",
    list4[[1]][2], " & ", list12[[1]][2], " \\\\ \n"
  ))

  cat(paste0(
    " & (", list5[[1]][2], ") & (", list13[[1]][2], ") & (",
    list6[[1]][2], ") & (", list14[[1]][2], ") & (",
    list7[[1]][2], ") & (", list15[[1]][2], ") & (",
    list8[[1]][2], ") & (", list16[[1]][2], ") \\\\ \n"
  ))

  cat(paste0(
    "$n=512$  & ", list1[[2]][3], " & ", list9[[2]][3], " & ",
    list2[[2]][3], " & ", list10[[2]][3], " & ",
    list3[[2]][3], " & ", list11[[2]][3], " & ",
    list4[[2]][3], " & ", list12[[2]][3], " \\\\ \n"
  ))

  cat(paste0(
    " & (", list5[[2]][3], ") & (", list13[[2]][3], ") & (",
    list6[[2]][3], ") & (", list14[[2]][3], ") & (",
    list7[[2]][3], ") & (", list15[[2]][3], ") & (",
    list8[[2]][3], ") & (", list16[[2]][3], ") \\\\ \n"
  ))

  cat(paste0(
    "$n=2048$  & ", list1[[3]][4], " & ", list9[[3]][4], " & ",
    list2[[3]][4], " & ", list10[[3]][4], " & ",
    list3[[3]][4], " & ", list11[[3]][4], " & ",
    list4[[3]][4], " & ", list12[[3]][4], " \\\\ \n"
  ))

  cat(paste0(
    " & (", list5[[3]][4], ") & (", list13[[3]][4], ") & (",
    list6[[3]][4], ") & (", list14[[3]][4], ") & (",
    list7[[3]][4], ") & (", list15[[3]][4], ") & (",
    list8[[3]][4], ") & (", list16[[3]][4], ") \\\\ \n"
  ))

  cat("\\end{tabular} \n")

  cat("\\end{table} \n")
}

latex_fun2(
  rejection_proportions_boot_A_1_indepground_depmethod,
  rejection_proportions_boot_B_1_indepground_depmethod,
  rejection_proportions_boot_C_1_indepground_depmethod,
  rejection_proportions_boot_D_1_indepground_depmethod,
  rejection_proportions_boot_a_1_5_indepground_depmethod,
  rejection_proportions_boot_b_1_5_indepground_depmethod,
  rejection_proportions_boot_c_1_5_indepground_depmethod,
  rejection_proportions_boot_d_1_5_indepground_depmethod,
  rejection_proportions_f_a_1_indepground,
  rejection_proportions_f_b_1_indepground,
  rejection_proportions_f_c_1_indepground,
  rejection_proportions_f_d_1_indepground,
  rejection_proportions_f_a_1_5_indepground,
  rejection_proportions_f_b_1_5_indepground,
  rejection_proportions_f_c_1_5_indepground,
  rejection_proportions_f_d_1_5_indepground
)


latex_fun2(
  rejection_proportions_boot_A_1_depground_depmethod,
  rejection_proportions_boot_B_1_depground_depmethod,
  rejection_proportions_boot_C_1_depground_depmethod,
  rejection_proportions_boot_D_1_depground_depmethod,
  rejection_proportions_boot_a_1_5_depground_depmethod,
  rejection_proportions_boot_b_1_5_depground_depmethod,
  rejection_proportions_boot_c_1_5_depground_depmethod,
  rejection_proportions_boot_d_1_5_depground_depmethod,
  rejection_proportions_f_a_1_depground,
  rejection_proportions_f_b_1_depground,
  rejection_proportions_f_c_1_depground,
  rejection_proportions_f_d_1_depground,
  rejection_proportions_f_a_1_5_depground,
  rejection_proportions_f_b_1_5_depground,
  rejection_proportions_f_c_1_5_depground,
  rejection_proportions_f_d_1_5_depground
)

################################################################################

# Illustrative plots

plot(model_a_sim(128), type = "l")


save.image(file.path(workspace_dir, "9_Comparison_sim_final.RData"))
