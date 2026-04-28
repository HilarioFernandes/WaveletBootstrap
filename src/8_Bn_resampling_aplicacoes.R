# =============================================================================
# 8_Bn_resampling_aplicacoes.R
# =============================================================================
# Purpose  : Real-data Bn applications (Nelson and Hohenthanner case studies).
# Chapter  : Chapter 4
# Inputs   : Nelson.txt, Hohenthanner.txt (Data files not included).
# Outputs  : Quasi U-statistic calculations and CI coverage for applications.
# Depends  : 7_Quasi_U_statistics_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
WORKSPACE_DIR <- file.path(BASE_PATH, "src", "WorkspaceData")
if (!dir.exists(WORKSPACE_DIR)) dir.create(WORKSPACE_DIR, recursive = TRUE)

source(file.path(BASE_PATH, "src", "7_Quasi_U_statistics_functions.R"))

# --- Testing Mode ---
TEST_MODE <- TRUE
# --------------------

# Define kernel used across all simulations
kernel <- function(x, y) (x - y)^2

B <- if (TEST_MODE) 5 else 100

################################################################################

#' Run Bn resampling tests (Bootstrap and Jackknife)
#'
#' @param X Data vector for group 1
#' @param Y Data vector for group 2
#' @param alpha Significance level
#' @param B Number of bootstrap reps
#' @return A list containing statistic/CI limits and p-values
Bn_resampling_tests <- function(X, Y, alpha, B) {
  n <- length(X) + length(Y)

  statistics_boot_discy <- rep(NA, B)

  statistics_boot_discn <- rep(NA, B)

  statistics_jackknife <- rep(NA, 2 * n)

  statistic <- B_n(X, Y, kernel)

  for (b in 1:B) {
    # Bootstrap with discrimination

    temp <- bootstrap_Bn(list(X, Y), TRUE)

    X_boot_discy <- temp[[1]]
    Y_boot_discy <- temp[[2]]

    statistics_boot_discy[b] <- B_n(X_boot_discy, Y_boot_discy, kernel)

    # Bootstrap without discrimination

    temp <- bootstrap_Bn(list(X, Y), FALSE)

    X_boot_discn <- temp[[1]]
    Y_boot_discn <- temp[[2]]

    statistics_boot_discn[b] <- B_n(X_boot_discn, Y_boot_discn, kernel)
  }

  # jackknife

  for (m in 1:(2 * n)) {
    if (m <= n) {
      statistics_jackknife[m] <- B_n(X[-m], Y, kernel)
    } else {
      statistics_jackknife[m] <- B_n(X, Y[-(m - n)], kernel)
    }
  }

  quantile_CI_boot_discy <- quantile(statistics_boot_discy, 1 - 0.05)
  quantile_CI_boot_discn <- quantile(statistics_boot_discn, 1 - 0.05)

  normal_approx_CI_boot_discy <- qnorm(1 - 0.05, mean = 0, sd = sd(statistics_boot_discy))

  normal_approx_CI_boot_discn <- qnorm(1 - 0.05, mean = 0, sd = sd(statistics_boot_discn))

  normal_approx_CI_jackknife <- qnorm(1 - 0.05,
    mean = 0,
    sd = sd(statistics_jackknife) * (2 * n - 1) / (sqrt(2 * n))
  )

  output1 <- c(
    statistic, quantile_CI_boot_discy, quantile_CI_boot_discn,
    normal_approx_CI_boot_discy, normal_approx_CI_boot_discn,
    normal_approx_CI_jackknife
  )

  names(output1) <- c(
    "statistic", "quant_boot_disc", "quant_boot_no_disc",
    "normal_boot_disc", "normal_boot_no_disc",
    "normal_jackknife"
  )

  output2 <- c(
    sum(statistics_boot_discy > statistic) / B,
    sum(statistics_boot_discn > statistic) / B,
    1 - pnorm(statistic, mean = 0, sd = sd(statistics_boot_discy)),
    1 - pnorm(statistic, mean = 0, sd = sd(statistics_boot_discn)),
    1 - pnorm(statistic,
      mean = 0,
      sd = sd(statistics_jackknife) * (2 * n - 1) / (sqrt(2 * n))
    )
  )

  return(list(output1, output2))
}

################################################################################

# Nelson

# Note: Data file is not included in the repository. See thesis for sources.
nelson <- read.csv(file.path(BASE_PATH, "Dados/Two sample tests", "Nelson.txt"), header = T)

plot(y = nelson[, 1], x = c(rep(1, 10), rep(2, 10)), xlim = c(0, 3))

X <- nelson[nelson[, 2] == 3, 1]

Y <- nelson[nelson[, 2] == 6, 1]

set.seed(0)

Bn_resampling_tests(X, Y, 0.05, 1000)

################################################################################

# hohenthanner

# Note: Data file is not included in the repository. See thesis for sources.
hohenthanner <- read.csv(file.path(BASE_PATH, "Dados/Two sample tests", "Hohenthanner.txt"), header = T)

plot(y = hohenthanner[, 1], x = c(rep(1, 5), rep(2, 5)), xlim = c(0, 3))

X <- hohenthanner[hohenthanner[, 2] == 1, 1]

Y <- hohenthanner[hohenthanner[, 2] == 2, 1]

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_aplicacoes_part_1.RData"))
set.seed(0)

Bn_resampling_tests(X, Y, 0.05, 1000)

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_aplicacoes_final.RData"))
