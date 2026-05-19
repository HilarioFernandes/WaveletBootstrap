# =============================================================================
# 5_CI_wv_aplicacoes.R
# =============================================================================
# Purpose  : Real-data CI applications for wavelet variance. (Ocean Shear, SPX/cboe)
# Chapter  : Chapter 3
# Inputs   : Ocean Shear.txt, SPX.csv (Data files not included in repository).
# Outputs  : Real-data CI plots (PNGs).
# Depends  : 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "2_Bootstrap_methods.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

library(multitaper)
b <- if (test_mode) 2 else 100

# Set and create output directory for plots
output_path <- file.path(base_path, "Plots/Plots_5")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

################################################################################

# this function calculates the wavelet variance standard deviation estimates
# z is the list with all squared wavelet coefficients (each list element corresponds to a scale)
# v is the list with all tapers (each list element corresponds to a scale)
sd_estimate_multitaper <- function(z, v) {
  n_levels <- length(z)

  n <- nrow(v[[1]])

  l <- (2^(1:n_levels) - 1) * (8 - 1) + 1

  # now we calculate the computable quantity J(0) and V(0)

  J <- matrix(0, ncol = n_levels, nrow = 5)
  V <- matrix(0, ncol = n_levels, nrow = 5)

  for (k in 1:5) {
    for (j in 1:n_levels) {
      # using non-boundary coeffs and the corresponding tapers
      # must set remove_boundary_coeffs <- TRUE
      J[k, j] <- z[[j]] %*% v[[j]][l[j]:n, k]


      # V[k,j] <- sum(v[[j]][l[j]:n,k])
      # V[k,j] <- sum(v[[j]][1:(n-l[j]+1),k])
      V[k, j] <- sum(v[[j]][, k])
    }
  }

  # calculating the estimated variance

  estimates <- NULL

  for (j in 1:n_levels) {
    # denominator of (12)

    denominator <- V[1, j]^2 + V[3, j]^2 + V[5, j]^2

    temp <- (J[1, j] * V[1, j] + J[3, j] * V[3, j] + J[5, j] * V[5, j]) / denominator

    sum <- 0

    for (k in 1:5) {
      sum <- sum + ((J[k, j] - temp * V[k, j])^2) / (n - l[j] + 1)
    }

    sum <- sum / 5

    estimates <- c(estimates, sum)
  }

  return(sqrt(estimates))
}

taper_fun <- function(n, l) {
  # return(dpss(n,5,3.5/(n-l+1))$v)

  # return(dpss(n,5,7/(n-l+1))$v)

  return(dpss(n, 5, 7)$v)
}


################################################################################

# Ocean Shear

# Note: Data file is not included in the repository. See thesis for sources.
ocean_shear <- read.table(file.path(base_path, "Dados", "Ocean Shear", "Ocean Shear.txt"))

ocean_shear <- cbind(1:nrow(ocean_shear), ocean_shear)

plot(ocean_shear, type = "l")

plot(diff(ocean_shear[, 2]), type = "l")

################################################################################

ci_wv_fun_global <- function(data, alpha, b) {
  n <- length(data)

  n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

  l <- (2^(1:n_levels) - 1) * (8 - 1) + 1

  wv_est <- wv_estimates(data)

  data_modwt <- modwt(data, n.levels = n_levels)

  wv_waveslim_nongaussian <- wave.variance(brick.wall(data_modwt, "la8", "modwt"),
    p = alpha / 2, type = "nongaussian"
  )[1:n_levels, ]

  wv_boot <- bootstrap_wavelet(data, "bw", TRUE, "SB", function(n) {
    1 / (4 * log2(n))
  }, b)

  # bootstrap variances estimates and point estimates for wavelet variances

  sd_boot <- apply(wv_boot, 2, function(x) {
    sd(x)
  })

  # assuming that the wavelet variance estimators are proportional to a chi-squared
  # distribution, we may estimate these two parameters via

  a_boot <- sd_boot^2 / (2 * wv_est)

  eta_boot <- (2 * wv_est^2) / sd_boot^2

  lower_boot <- a_boot * qchisq(rep(alpha / 2, n_levels), df = eta_boot)
  upper_boot <- a_boot * qchisq(rep(1 - (alpha / 2), n_levels), df = eta_boot)

  # multitaper confidence intervals (manually implemented)

  taper_fun_temp <- function(l) {
    return(taper_fun(n, l))
  }

  v <- lapply(l, taper_fun_temp)

  sd_multitaper <- sd_estimate_multitaper(nonbnd_sq_modwt_coeffs(data), v)

  a_multitaper <- sd_multitaper^2 / (2 * wv_est)
  eta_multitaper <- (2 * wv_est^2) / sd_multitaper^2

  lower_multitaper <- a_multitaper * qchisq(rep(alpha / 2, n_levels), df = eta_multitaper)
  upper_multitaper <- a_multitaper * qchisq(rep(1 - (alpha / 2), n_levels), df = eta_multitaper)

  # output

  output <- cbind(wv_waveslim_nongaussian, wv_est, lower_boot, upper_boot, lower_multitaper, upper_multitaper)

  colnames(output) <- c(
    "wv_periodic",
    "lower_multitaper",
    "upper_multitaper",
    "wv_nonbound",
    "lower_boot",
    "upper_boot",
    "lower_multitaper_2",
    "upper_multitaper_2"
  )

  return(output)
}

ci_wv_fun_global(ocean_shear[, 2], 0.05, 100)

ci_wv_fun_local <- function(data, windowsize, alpha, b) {
  n <- length(data)

  output <- vector(mode = "list", length = n - windowsize + 1)

  for (i in 1:(n - windowsize + 1)) {
    message(round(i / (n - windowsize + 1), 3), "\r", appendLF = FALSE)
    flush.console()
    Sys.sleep(0.1)

    data_window <- data[i:(i + windowsize - 1)]

    output[[i]] <- ci_wv_fun_global(data_window, alpha, b)
  }

  return(output)
}

set.seed(0)
results_ocean_shear <- ci_wv_fun_local(diff(ocean_shear[150:2250, 2]), 512, 0.05, b)

col1 <- "grey"
col2 <- "black"

lty1 <- "solid"
lty2 <- "dotted"

{
  png(
    file = file.path(output_path, "5_OceanShear.png"),
    width = 1500, height = 1000, res = 150
  )

  par(mfrow = c(2, 2))
  par(mar = c(4.1, 4.1, 1.5, 0.5))

  {
    plot(
      unlist(lapply(results_ocean_shear, function(x) {
        x[1, 1]
      })),
      type = "l", ylim = c(0, 0.0016),
      ylab = "var. de ondaletas (j=1)", xaxt = "n", yaxt = "n", xlab = "",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.0016, length.out = 3), labels = c("0", "", "1.6e-3"))
    axis(side = 1, at = seq(1, 1589, length.out = 7), labels = F)

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[1, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[1, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[1, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[1, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )

    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[1,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[1,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_ocean_shear, function(x) {
        x[2, 1]
      })),
      type = "l", ylim = c(0, 0.0024),
      ylab = "var. de ondaletas (j=2)", xaxt = "n", yaxt = "n", xlab = "",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.0024, length.out = 3), labels = c("0", "", "2.4e-3"))
    axis(side = 1, at = seq(1, 1589, length.out = 7), labels = F)

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[2, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[2, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[2, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[2, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[2,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[2,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_ocean_shear, function(x) {
        x[3, 1]
      })),
      type = "l", ylim = c(0, 0.0018),
      ylab = "var. de ondaletas (j=3)", xaxt = "n", yaxt = "n", xlab = "Profundidade (metros)",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.0018, length.out = 3), labels = c("0", "", "1.8e-3"))
    axis(side = 1, at = seq(1, 1589, length.out = 7), labels = seq(365, 575, length.out = 7))

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[3, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[3, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[3, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[3, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[3,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[3,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_ocean_shear, function(x) {
        x[4, 1]
      })),
      type = "l", ylim = c(0, 0.0012),
      ylab = "var. de ondaletas (j=4)", xaxt = "n", yaxt = "n", xlab = "Profundidade (metros)",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.0012, length.out = 3), labels = c("0", "", "1.2e-3"))
    axis(side = 1, at = seq(1, 1589, length.out = 7), labels = seq(365, 575, length.out = 7))

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[4, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[4, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[4, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_ocean_shear, function(x) {
        x[4, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[4,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[4,3]})), lty = "dashed")
  }

  dev.off()
}

################################################################################

# cboe

# Function to calculate returns
returns_fun <- function(data) {
  data_temp <- unname(unlist(as.vector(data[, 2:ncol(data)])))

  output <- NULL

  for (i in 2:length(data_temp)) {
    output <- c(output, (data_temp[i] - data_temp[i - 1]) / data_temp[i - 1])
  }

  return(output)
}

# Function to subset the data. It returns a matrix of indexes.
# interval_type is either "minutes", "hours" or "days".
# interval is the amount of interval_type between samples (except if interval_type == "days",
# in this case interval is ignored).
subset_data <- function(interval_type, interval, day_start, day_stop) {
  if (interval_type == "minutes") {
    adj_interval <- 60 * interval
  } else if (interval_type == "hours") {
    adj_interval <- 3600 * interval
  } else if (interval_type == "days") {
    adj_interval <- 24600
  }

  output1 <- seq(1, 24600, by = adj_interval)

  return(list(output1, c(1, day_start:day_stop)))
}


# Note: Data file is not included in the repository. See thesis for sources.
cboe <- read.csv(file.path(base_path, "Dados", "SPX_second", "SPX.csv"), check.names = FALSE)

# sampling
indexes <- subset_data("minutes", 5, 12, 26)

cboe_selection <- cboe[indexes[[1]], indexes[[2]]]

cboe_selection_returns <- returns_fun(cboe_selection)

plot(cboe_selection_returns, type = "l")

save.image(file.path(workspace_dir, "5_CI_wv_aplicacoes_part_1.RData"))
set.seed(0)

results_cboe <- ci_wv_fun_local(cboe_selection_returns, 512, 0.05, b)

labels_plot <- rep(substr(colnames(cboe_selection)[seq(2, ncol(cboe_selection), by = 1)], 9, 10), each = 82)
# labels_plot <- as.integer(labels_plot[513:1230])
labels_plot <- as.integer(labels_plot[1:718])

min_indexes_labels <- NULL

for (i in unique(labels_plot)) {
  min_indexes_labels <- c(min_indexes_labels, which(labels_plot == i)[1])
}

{
  png(
    file = file.path(output_path, "5_CBOE.png"),
    width = 1500, height = 1000, res = 150
  )

  par(mfrow = c(2, 2))
  par(mar = c(4.1, 4.1, 1.5, 0.5))

  {
    plot(
      unlist(lapply(results_cboe, function(x) {
        x[1, 1]
      })),
      type = "l", ylim = c(0, 0.00006),
      ylab = "var. de ondaletas (j=1)", xaxt = "n", yaxt = "n", xlab = "",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.00006, length.out = 3), labels = c("0", "", "6e-5"))
    axis(
      side = 1, at = min_indexes_labels,
      labels = F
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[1, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[1, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[1, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[1, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[1,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[1,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_cboe, function(x) {
        x[2, 1]
      })),
      type = "l", ylim = c(0, 0.00003),
      ylab = "var. de ondaletas (j=2)", yaxt = "n", xaxt = "n", xlab = "",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.00003, length.out = 3), labels = c("0", "", "3e-5"))
    axis(
      side = 1, at = min_indexes_labels,
      labels = F
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[2, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[2, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[2, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[2, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[2,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[2,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_cboe, function(x) {
        x[3, 1]
      })),
      type = "l", ylim = c(0, 0.00002),
      ylab = "var. de ondaletas (j=3)", xaxt = "n", yaxt = "n", xlab = "Dia",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.00002, length.out = 3), labels = c("0", "", "2e-5"))
    axis(
      side = 1, at = min_indexes_labels,
      labels = labels_plot[min_indexes_labels]
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[3, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[3, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[3, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[3, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[3,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[3,3]})), lty = "dashed")
  }

  {
    plot(
      unlist(lapply(results_cboe, function(x) {
        x[4, 1]
      })),
      type = "l", ylim = c(0, 0.000012),
      ylab = "var. de ondaletas (j=4)", xaxt = "n", yaxt = "n", xlab = "Dia",
      cex.lab = 1.25, cex.main = 1.25, cex.axis = 1.25, cex.sub = 1.25,
      lwd = 1.5
    )
    axis(side = 2, at = seq(0, 0.000012, length.out = 3), labels = c("0", "", "1.2e-5"))
    axis(
      side = 1, at = min_indexes_labels,
      labels = labels_plot[min_indexes_labels]
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[4, 5]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[4, 6]
      })),
      col = col1, lty = lty1,
      lwd = 1.5
    )

    lines(
      unlist(lapply(results_cboe, function(x) {
        x[4, 7]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    lines(
      unlist(lapply(results_cboe, function(x) {
        x[4, 8]
      })),
      col = col2, lty = lty2,
      lwd = 1.5
    )
    legend("topright",
      legend = c("Est. local", "Bootstrap", "Multitaper"),
      col = c("black", col1, col2),
      lty = c("solid", lty1, lty2),
      lwd = 1.5
    )

    # lines(unlist(lapply(results_cboe, function(x){x[4,2]})), lty = "dashed")
    # lines(unlist(lapply(results_cboe, function(x){x[4,3]})), lty = "dashed")
  }
  dev.off()
}

save.image(file.path(workspace_dir, "5_CI_wv_aplicacoes_final.RData"))
