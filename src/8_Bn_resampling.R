# =============================================================================
# 8_Bn_resampling.R
# =============================================================================
# Purpose  : Simulation study of Bn distribution via resampling methods.
# Chapter  : Chapter 4, Appendix C.1
# Inputs   : 7_Quasi_U_statistics_functions.R
# Outputs  : Size and power study results for the Bn statistic.
# Depends  : 7_Quasi_U_statistics_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "7_Quasi_U_statistics_functions.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

# Define kernel used across all simulations
kernel <- function(x, y) (x - y)^2

n <- 100
iterations <- if (test_mode) 2 else 100
b <- if (test_mode) 5 else 100


################################################################################

x_e <- vector(mode = "list", length = iterations)
y_e <- vector(mode = "list", length = iterations)

statistics_e <- NULL

statistics_boot_discy_e <- matrix(NA, nrow = iterations, ncol = b)

statistics_boot_discn_e <- matrix(NA, nrow = iterations, ncol = b)

statistics_jackknife_e <- matrix(NA, nrow = iterations, ncol = 2 * n)


set.seed(46)

for (i in 1:iterations) {
  x_e[[i]] <- rnorm(n)
  y_e[[i]] <- rnorm(n)

  statistics_e <- c(statistics_e, b_n(x_e[[i]], y_e[[i]], kernel))
}

set.seed(47)

for (i in 1:iterations) {
  x <- x_e[[i]]
  y <- y_e[[i]]

  for (b in 1:b) {
    # Bootstrap with discrimination

    temp <- bootstrap_bn(list(x, y), TRUE)

    x_boot_discy <- temp[[1]]
    y_boot_discy <- temp[[2]]

    statistics_boot_discy_e[i, b] <- b_n(x_boot_discy, y_boot_discy, kernel)

    # Bootstrap without discrimination

    temp <- bootstrap_bn(list(x, y), FALSE)

    x_boot_discn <- temp[[1]]
    y_boot_discn <- temp[[2]]

    statistics_boot_discn_e[i, b] <- b_n(x_boot_discn, y_boot_discn, kernel)
  }
}


set.seed(48)

for (i in 1:iterations) {
  x <- x_e[[i]]
  y <- y_e[[i]]

  # jackknife

  for (m in 1:(2 * n)) {
    if (m <= n) {
      statistics_jackknife_e[i, m] <- b_n(x[-m], y, kernel)
    } else {
      statistics_jackknife_e[i, m] <- b_n(x, y[-(m - n)], kernel)
    }
  }
}

save.image(file.path(workspace_dir, "8_Bn_resampling_e.RData"))

################################################################################

x_f <- vector(mode = "list", length = iterations)
y_f <- vector(mode = "list", length = iterations)

statistics_f <- NULL

statistics_boot_discy_f <- matrix(NA, nrow = iterations, ncol = b)

statistics_boot_discn_f <- matrix(NA, nrow = iterations, ncol = b)

statistics_jackknife_f <- matrix(NA, nrow = iterations, ncol = 2 * n)


set.seed(49)

for (i in 1:iterations) {
  x_f[[i]] <- rnorm(n)
  y_f[[i]] <- rnorm(n, 1, 1)

  statistics_f <- c(statistics_f, b_n(x_f[[i]], y_f[[i]], kernel))
}


set.seed(50)

for (i in 1:iterations) {
  x <- x_f[[i]]
  y <- y_f[[i]]

  for (b in 1:b) {
    # Bootstrap with discrimination

    temp <- bootstrap_bn(list(x, y), TRUE)

    x_boot_discy <- temp[[1]]
    y_boot_discy <- temp[[2]]

    statistics_boot_discy_f[i, b] <- b_n(x_boot_discy, y_boot_discy, kernel)

    # Bootstrap without discrimination

    temp <- bootstrap_bn(list(x, y), FALSE)

    x_boot_discn <- temp[[1]]
    y_boot_discn <- temp[[2]]

    statistics_boot_discn_f[i, b] <- b_n(x_boot_discn, y_boot_discn, kernel)
  }
}


set.seed(51)

for (i in 1:iterations) {
  x <- x_f[[i]]
  y <- y_f[[i]]

  # jackknife

  for (m in 1:(2 * n)) {
    if (m <= n) {
      statistics_jackknife_f[i, m] <- b_n(x[-m], y, kernel)
    } else {
      statistics_jackknife_f[i, m] <- b_n(x, y[-(m - n)], kernel)
    }
  }
}

save.image(file.path(workspace_dir, "8_Bn_resampling_f.RData"))

################################################################################

x_g <- vector(mode = "list", length = iterations)
y_g <- vector(mode = "list", length = iterations)

statistics_g <- NULL

statistics_boot_discy_g <- matrix(NA, nrow = iterations, ncol = b)

statistics_boot_discn_g <- matrix(NA, nrow = iterations, ncol = b)

statistics_jackknife_g <- matrix(NA, nrow = iterations, ncol = 2 * n)


set.seed(52)

for (i in 1:iterations) {
  x_g[[i]] <- rchisq(n, df = 3)
  y_g[[i]] <- rchisq(n, df = 3)

  statistics_g <- c(statistics_g, b_n(x_g[[i]], y_g[[i]], kernel))
}


set.seed(53)

for (i in 1:iterations) {
  x <- x_g[[i]]
  y <- y_g[[i]]

  for (b in 1:b) {
    # Bootstrap with discrimination

    temp <- bootstrap_bn(list(x, y), TRUE)

    x_boot_discy <- temp[[1]]
    y_boot_discy <- temp[[2]]

    statistics_boot_discy_g[i, b] <- b_n(x_boot_discy, y_boot_discy, kernel)

    # Bootstrap without discrimination

    temp <- bootstrap_bn(list(x, y), FALSE)

    x_boot_discn <- temp[[1]]
    y_boot_discn <- temp[[2]]

    statistics_boot_discn_g[i, b] <- b_n(x_boot_discn, y_boot_discn, kernel)
  }
}

set.seed(54)

for (i in 1:iterations) {
  x <- x_g[[i]]
  y <- y_g[[i]]

  # jackknife

  for (m in 1:(2 * n)) {
    if (m <= n) {
      statistics_jackknife_g[i, m] <- b_n(x[-m], y, kernel)
    } else {
      statistics_jackknife_g[i, m] <- b_n(x, y[-(m - n)], kernel)
    }
  }
}

save.image(file.path(workspace_dir, "8_Bn_resampling_g.RData"))

################################################################################


x_h <- vector(mode = "list", length = iterations)
y_h <- vector(mode = "list", length = iterations)

statistics_h <- NULL

statistics_boot_discy_h <- matrix(NA, nrow = iterations, ncol = b)

statistics_boot_discn_h <- matrix(NA, nrow = iterations, ncol = b)

statistics_jackknife_h <- matrix(NA, nrow = iterations, ncol = 2 * n)


set.seed(55)

for (i in 1:iterations) {
  x_h[[i]] <- rchisq(n, df = 3)
  y_h[[i]] <- rchisq(n, df = 6)

  statistics_h <- c(statistics_h, b_n(x_h[[i]], y_h[[i]], kernel))
}


set.seed(56)

for (i in 1:iterations) {
  x <- x_h[[i]]
  y <- y_h[[i]]

  for (b in 1:b) {
    # Bootstrap with discrimination

    temp <- bootstrap_bn(list(x, y), TRUE)

    x_boot_discy <- temp[[1]]
    y_boot_discy <- temp[[2]]

    statistics_boot_discy_h[i, b] <- b_n(x_boot_discy, y_boot_discy, kernel)

    # Bootstrap without discrimination

    temp <- bootstrap_bn(list(x, y), FALSE)

    x_boot_discn <- temp[[1]]
    y_boot_discn <- temp[[2]]

    statistics_boot_discn_h[i, b] <- b_n(x_boot_discn, y_boot_discn, kernel)
  }
}

set.seed(57)

for (i in 1:iterations) {
  x <- x_h[[i]]
  y <- y_h[[i]]

  # jackknife

  for (m in 1:(2 * n)) {
    if (m <= n) {
      statistics_jackknife_h[i, m] <- b_n(x[-m], y, kernel)
    } else {
      statistics_jackknife_h[i, m] <- b_n(x, y[-(m - n)], kernel)
    }
  }
}


save.image(file.path(workspace_dir, "8_Bn_resampling_h.RData"))

################################################################################

quantile_ci_boot_discy_e <- quantile_ci(statistics_boot_discy_e, 0.05)
quantile_ci_boot_discn_e <- quantile_ci(statistics_boot_discn_e, 0.05)

normal_approx_ci_boot_discy_e <- approx_ci(
  statistics_boot_discy_e,
  "boot", 0.05
)

normal_approx_ci_boot_discn_e <- approx_ci(
  statistics_boot_discn_e,
  "boot", 0.05
)

normal_approx_ci_jackknife_e <- approx_ci(
  statistics_jackknife_e,
  "jackknife", 0.05
)


quantile_ci_boot_discy_f <- quantile_ci(statistics_boot_discy_f, 0.05)
quantile_ci_boot_discn_f <- quantile_ci(statistics_boot_discn_f, 0.05)

normal_approx_ci_boot_discy_f <- approx_ci(
  statistics_boot_discy_f,
  "boot", 0.05
)

normal_approx_ci_boot_discn_f <- approx_ci(
  statistics_boot_discn_f,
  "boot", 0.05
)

normal_approx_ci_jackknife_f <- approx_ci(
  statistics_jackknife_f,
  "jackknife", 0.05
)


quantile_ci_boot_discy_g <- quantile_ci(statistics_boot_discy_g, 0.05)
quantile_ci_boot_discn_g <- quantile_ci(statistics_boot_discn_g, 0.05)

normal_approx_ci_boot_discy_g <- approx_ci(
  statistics_boot_discy_g,
  "boot", 0.05
)

normal_approx_ci_boot_discn_g <- approx_ci(
  statistics_boot_discn_g,
  "boot", 0.05
)

normal_approx_ci_jackknife_g <- approx_ci(
  statistics_jackknife_g,
  "jackknife", 0.05
)


quantile_ci_boot_discy_h <- quantile_ci(statistics_boot_discy_h, 0.05)
quantile_ci_boot_discn_h <- quantile_ci(statistics_boot_discn_h, 0.05)

normal_approx_ci_boot_discy_h <- approx_ci(
  statistics_boot_discy_h,
  "boot", 0.05
)

normal_approx_ci_boot_discn_h <- approx_ci(
  statistics_boot_discn_h,
  "boot", 0.05
)

normal_approx_ci_jackknife_h <- approx_ci(
  statistics_jackknife_h,
  "jackknife", 0.05
)


################################################################################


table_results <- matrix(NA, nrow = 4, ncol = 5)
colnames(table_results) <- c(
  "quantile_boot_disc",
  "quantile_boot",
  "normal_boot_disc",
  "normal_boot",
  "normal_jack"
)
rownames(table_results) <- c("e", "f", "g", "h")


table_results[1, ] <- c(
  coverage_ci(statistics_e, t(quantile_ci_boot_discy_e)),
  coverage_ci(statistics_e, t(quantile_ci_boot_discn_e)),
  coverage_ci(statistics_e, t(normal_approx_ci_boot_discy_e)),
  coverage_ci(statistics_e, t(normal_approx_ci_boot_discn_e)),
  coverage_ci(statistics_e, t(normal_approx_ci_jackknife_e))
)

table_results[2, ] <- c(
  coverage_ci(statistics_f, t(quantile_ci_boot_discy_f)),
  coverage_ci(statistics_f, t(quantile_ci_boot_discn_f)),
  coverage_ci(statistics_f, t(normal_approx_ci_boot_discy_f)),
  coverage_ci(statistics_f, t(normal_approx_ci_boot_discn_f)),
  coverage_ci(statistics_f, t(normal_approx_ci_jackknife_f))
)

table_results[3, ] <- c(
  coverage_ci(statistics_g, t(quantile_ci_boot_discy_g)),
  coverage_ci(statistics_g, t(quantile_ci_boot_discn_g)),
  coverage_ci(statistics_g, t(normal_approx_ci_boot_discy_g)),
  coverage_ci(statistics_g, t(normal_approx_ci_boot_discn_g)),
  coverage_ci(statistics_g, t(normal_approx_ci_jackknife_g))
)

table_results[4, ] <- c(
  coverage_ci(statistics_h, t(quantile_ci_boot_discy_h)),
  coverage_ci(statistics_h, t(quantile_ci_boot_discn_h)),
  coverage_ci(statistics_h, t(normal_approx_ci_boot_discy_h)),
  coverage_ci(statistics_h, t(normal_approx_ci_boot_discn_h)),
  coverage_ci(statistics_h, t(normal_approx_ci_jackknife_h))
)

View(table_results)

save.image(file.path(workspace_dir, "8_Bn_resampling_final.RData"))
