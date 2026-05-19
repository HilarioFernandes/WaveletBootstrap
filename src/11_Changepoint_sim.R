# =============================================================================
# 11_Changepoint_sim.R
# =============================================================================
# Purpose  : Simulation study for changepoint detection performance.
# Chapter  : Chapter 3, Appendix b.5
# Inputs   : 10_Comparison_Changepoint_functions.R
# Outputs  : Power and size rejection rates for changepoint detection.
# Depends  : 10_Comparison_Changepoint_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "10_Comparison_Changepoint_functions.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

################################################################################

#' Run changepoint detection simulation (n=2048, block=128)
#'
#' @param model Model identifier ("A" through "D")
#' @param iterations Number of simulation iterations
#' @param b Number of bootstrap reps
#' @param alpha_0 Target significance level
#' @param change Logical. If TRUE, introduces a variance change.
#' @param change_type Either "abrupt" or "drift"
#' @return Matrix of rejection proportions per block
changepoint_sim <- function(model, iterations, b, alpha_0, change, change_type) {
  sim_fun <- get(paste0("model_", model, "_sim"))

  results <- matrix(0, nrow = 15, ncol = 4)

  rownames(results) <- 1:15

  for (i in 1:iterations) {
    print(i)

    x <- sim_fun(2048)

    if (change == TRUE) {
      if (change_type == "abrupt") {
        x[501:1471] <- 1.25 * x[501:1471]
        x[1472:2048] <- (1.25^2) * x[1472:2048]
      } else if (change_type == "drift") {
        x <- seq(1, 1.25^2, length.out = length(x)) * x
      }
    }

    # plot(x, type = "l")

    temp <- changepointdetection(x, 128, FALSE, "F", alpha_0, 2, NA)
    results[temp[, 2] - 1, 1] <- results[temp[, 2] - 1, 1] + 1

    temp <- changepointdetection(x, 128, TRUE, "F", alpha_0, 2, NA)
    results[temp[, 2] - 1, 2] <- results[temp[, 2] - 1, 2] + 1

    temp <- changepointdetection(x, 128, FALSE, "boot_quant", alpha_0, 2, b)
    results[temp[, 2] - 1, 3] <- results[temp[, 2] - 1, 3] + 1

    temp <- changepointdetection(x, 128, TRUE, "boot_quant", alpha_0, 2, b)
    results[temp[, 2] - 1, 4] <- results[temp[, 2] - 1, 4] + 1
  }

  return(results / iterations)
}

changepoint_sim_b_abrupt <- changepoint_sim("B", 10, 100, 0.05, TRUE, "abrupt")
changepoint_sim_b_drift <- changepoint_sim("B", 10, 100, 0.05, TRUE, "drift")
changepoint_sim_b_na <- changepoint_sim("B", 10, 100, 0.05, FALSE, NA)

################################################################################

iterations <- if (test_mode) 2 else 100
b <- if (test_mode) 2 else 100

################################################################################
# 128

set.seed(66)

start_time <- Sys.time()

changepoint_sim_a_abrupt <- changepoint_sim("A", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_a_drift <- changepoint_sim("A", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_a_na <- changepoint_sim("A", iterations, b, 0.05, FALSE, NA)

changepoint_sim_b_abrupt <- changepoint_sim("B", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_b_drift <- changepoint_sim("B", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_b_na <- changepoint_sim("B", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "11_Changepoint_sim_part_1.RData"))
set.seed(67)

start_time <- Sys.time()

changepoint_sim_c_abrupt <- changepoint_sim("C", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_c_drift <- changepoint_sim("C", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_c_na <- changepoint_sim("C", iterations, b, 0.05, FALSE, NA)

changepoint_sim_d_abrupt <- changepoint_sim("D", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_d_drift <- changepoint_sim("D", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_d_na <- changepoint_sim("D", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

################################################################################

#' Run changepoint detection simulation (n=8192, block=512)
#'
#' @param model Model identifier ("A" through "D")
#' @param iterations Number of simulation iterations
#' @param b Number of bootstrap reps
#' @param alpha_0 Target significance level
#' @param change Logical. If TRUE, introduces a variance change.
#' @param change_type Either "abrupt" or "drift"
#' @return Matrix of rejection proportions per block
changepoint_sim_512 <- function(model, iterations, b, alpha_0, change, change_type) {
  sim_fun <- get(paste0("model_", model, "_sim"))

  results <- matrix(0, nrow = 15, ncol = 4)

  rownames(results) <- 1:15

  for (i in 1:iterations) {
    print(i)

    x <- sim_fun(8192)

    if (change == TRUE) {
      if (change_type == "abrupt") {
        x[2004:5884] <- 1.25 * x[2004:5884]
        x[5885:8192] <- (1.25^2) * x[5885:8192]
      } else if (change_type == "drift") {
        x <- seq(1, 1.25^2, length.out = length(x)) * x
      }
    }

    # plot(x, type = "l")

    temp <- changepointdetection(x, 512, FALSE, "F", alpha_0, 3, NA)
    results[temp[, 2] - 1, 1] <- results[temp[, 2] - 1, 1] + 1

    temp <- changepointdetection(x, 512, TRUE, "F", alpha_0, 3, NA)
    results[temp[, 2] - 1, 2] <- results[temp[, 2] - 1, 2] + 1

    temp <- changepointdetection(x, 512, FALSE, "boot_quant", alpha_0, 3, b)
    results[temp[, 2] - 1, 3] <- results[temp[, 2] - 1, 3] + 1

    temp <- changepointdetection(x, 512, TRUE, "boot_quant", alpha_0, 3, b)
    results[temp[, 2] - 1, 4] <- results[temp[, 2] - 1, 4] + 1
  }

  return(results / iterations)
}

changepoint_sim_512("B", 10, 100, 0.05, TRUE, "abrupt")

################################################################################
# 512

save.image(file.path(workspace_dir, "11_Changepoint_sim_part_2.RData"))
set.seed(68)

start_time <- Sys.time()

changepoint_sim_a_abrupt_512 <- changepoint_sim_512("A", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_a_drift_512 <- changepoint_sim_512("A", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_a_na_512 <- changepoint_sim_512("A", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "11_Changepoint_sim_part_3.RData"))
set.seed(69)

start_time <- Sys.time()

changepoint_sim_b_abrupt_512 <- changepoint_sim_512("B", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_b_drift_512 <- changepoint_sim_512("B", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_b_na_512 <- changepoint_sim_512("B", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
  print(duration)

save.image(file.path(workspace_dir, "11_Changepoint_sim_part_4.RData"))
set.seed(70)

start_time <- Sys.time()

changepoint_sim_c_abrupt_512 <- changepoint_sim_512("C", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_c_drift_512 <- changepoint_sim_512("C", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_c_na_512 <- changepoint_sim_512("C", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

save.image(file.path(workspace_dir, "11_Changepoint_sim_part_5.RData"))
set.seed(71)

start_time <- Sys.time()

changepoint_sim_d_abrupt_512 <- changepoint_sim_512("D", iterations, b, 0.05, TRUE, "abrupt")
changepoint_sim_d_drift_512 <- changepoint_sim_512("D", iterations, b, 0.05, TRUE, "drift")
changepoint_sim_d_na_512 <- changepoint_sim_512("D", iterations, b, 0.05, FALSE, NA)

end_time <- Sys.time()

duration <- end_time - start_time
print(duration)

################################################################################

#' Generate LaTeX table for changepoint simulation results
#'
#' @param df Results matrix for first sample size (e.g., 128)
#' @param df2 Results matrix for second sample size (e.g., 512)
#' @return Prints a formatted LaTeX table to the console
latex_fun <- function(df, df2) {
  cat("\\begin{table}[h!] \n")

  cat("\\centering \n")

  cat("\\caption{Taxas empíricas de rejeição para $l=$ e modelo () com/sem mudança abrupta/gradual de variância.} \n")

  cat("\\label{temp} \n")

  cat("\\begin{tabular}{c|cccc} \n")

  cat("Update & Aprox. $F$ & Aprox. $F$ (seq.) & Bootstrap & Bootstrap (seq.) \\\\ \\hline \n")

  for (i in 1:15) {
    cat(paste0(i, " & ", paste(df[i, ], collapse = " & "), " \\\\ \n"))
    cat(paste0(" & (", paste(df2[i, ], collapse = ") & ("), ") \\\\ \\hline \n"))
  }

  cat("\\end{tabular} \n")

  cat("\\end{table} \n")
}

latex_fun(changepoint_sim_a_na, changepoint_sim_a_na_512)
latex_fun(changepoint_sim_a_drift, changepoint_sim_a_drift_512)
latex_fun(changepoint_sim_a_abrupt, changepoint_sim_a_abrupt_512)

latex_fun(changepoint_sim_b_na, changepoint_sim_b_na_512)
latex_fun(changepoint_sim_b_drift, changepoint_sim_b_drift_512)
latex_fun(changepoint_sim_b_abrupt, changepoint_sim_b_abrupt_512)

latex_fun(changepoint_sim_c_na, changepoint_sim_c_na_512)
latex_fun(changepoint_sim_c_drift, changepoint_sim_c_drift_512)
latex_fun(changepoint_sim_c_abrupt, changepoint_sim_c_abrupt_512)

latex_fun(changepoint_sim_d_na, changepoint_sim_d_na_512)
latex_fun(changepoint_sim_d_drift, changepoint_sim_d_drift_512)
latex_fun(changepoint_sim_d_abrupt, changepoint_sim_d_abrupt_512)


save.image(file.path(workspace_dir, "11_Changepoint_sim_final.RData"))
