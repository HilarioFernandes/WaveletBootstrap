# =============================================================================
# 11_Changepoint_sim_aplicacoes.R
# =============================================================================
# Purpose  : Real-data changepoint applications (Ocean Shear, cboe SPX).
# Chapter  : Chapter 3
# Inputs   : Ocean Shear.txt, SPX.csv (Data files not included).
# Outputs  : Detected changepoints and visualization of variance shifts.
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

library(dplyr)

# Set and create output directory for plots
output_path <- file.path(base_path, "Plots/Plots_11")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


b <- if (test_mode) 5 else 1000
alpha <- 0.05

################################################################################

#' Run all changepoint detection variants
#'
#' @param data Time series data vector
#' @param block_size Length of each block
#' @param alpha_0 Significance level
#' @param jmax Maximum wavelet level
#' @param b Number of bootstrap reps
#' @return A list with detection results for (F, F-seq, Boot, Boot-seq)
changepointdetection_fun <- function(data, block_size, alpha_0, jmax, b) {
  output <- list(NULL, NULL, NULL, NULL)

  output[[1]] <- changepointdetection(data, block_length, FALSE, "F", alpha_0, jmax, b)
  output[[2]] <- changepointdetection(data, block_length, TRUE, "F", alpha_0, jmax, b)
  output[[3]] <- changepointdetection(data, block_length, FALSE, "boot_quant", alpha_0, jmax, b)
  output[[4]] <- changepointdetection(data, block_length, TRUE, "boot_quant", alpha_0, jmax, b)

  return(output)
}

################################################################################

# Note: Data file is not included in the repository. See thesis for sources.
ocean_shear <- unname(unlist(as.vector(read.table(file.path(base_path, "Dados", "Ocean Shear", "Ocean Shear.txt")))))

plot(ocean_shear, type = "l")

x_oceanshear <- diff(ocean_shear)[1:(512 * 13)]

block_length <- 512
jmax <- 3

meter_indexes <- seq(350.1, by = 0.1, length.out = length(x_oceanshear))

{
  png(file = file.path(output_path, "OceanShearApplication.png"), width = 1800, height = 900, res = 210)

  par(mar = c(4, 5, 2, 2), mfrow = c(1, 1))
  plot(x_oceanshear,
    type = "l", xaxt = "n", xlab = "Profundidade (metros)",
    ylab = "Cisalhamento \n (diferenciado)",
    lwd = 1.5
  )
  abline(
    v = seq(0.5, 6656.5, by = block_length), lty = "dashed",
    lwd = 1.5
  )
  axis(1,
    at = seq(1, length(x_oceanshear), length.out = 14),
    labels = round(seq(350.1, 350 + length(x_oceanshear) / 10, length.out = 14), 0)
  )

  dev.off()
}

set.seed(0)
results_oceanshear <- changepointdetection_fun(x_oceanshear, block_length, alpha, jmax, 100)

################################################################################

#' Calculate relative returns from price data
#'
#' @param data Dataframe with prices (first column dates/times)
#' @return Vector of relative returns
returns_fun <- function(data) {
  data_temp <- unname(unlist(as.vector(data[, 2:ncol(data)])))

  output <- (data_temp[2:length(data_temp)] - data_temp[1:(length(data_temp) - 1)]) / data_temp[1:(length(data_temp) - 1)]

  names_temp <- expand.grid(data[, 1], colnames(data)[2:ncol(data)])


  names_temp <- apply(names_temp, 1, function(x) {
    paste0(x[2], "_", x[1])
  })

  names(output) <- names_temp[2:length(names_temp)]

  return(output)
}

#' Subset data by time intervals
#'
#' @param interval_type Either "minutes", "hours", "days", or "seconds"
#' @param interval Step size for subsetting
#' @param day_start Start index for days
#' @param day_stop Stop index for days
#' @return A list with time indexes and day selection indexes
subset_data <- function(interval_type, interval, day_start, day_stop) {
  if (interval_type == "minutes") {
    adj_interval <- 60 * interval
  } else if (interval_type == "hours") {
    adj_interval <- 3600 * interval
  } else if (interval_type == "days") {
    adj_interval <- 24600
  } else {
    adj_interval <- interval
  }

  output1 <- seq(1, 24600, by = adj_interval)

  return(list(output1, c(1, day_start:day_stop)))
}

# Note: Data file is not included in the repository. See thesis for sources.
cboe <- read.csv(file.path(base_path, "Dados", "SPX_second", "SPX.csv"), check.names = FALSE)

# sampling
indexes <- subset_data("seconds", 5, 2, 2)

cboe_selection <- cboe[indexes[[1]], indexes[[2]]]

cboe_selection_returns <- returns_fun(cboe_selection)

x_cboe <- cboe_selection_returns[1:2048]

block_length <- 128
jmax <- 2

{
  png(file = file.path(output_path, "CBOEApplication.png"), width = 1800, height = 900, res = 210)
  par(mar = c(4, 4, 2, 2))
  plot(x_cboe,
    type = "l", xaxt = "n", xlab = "Horário (hh:mm)", ylab = "Retorno",
    ylim = c(-0.001, 0.001),
    lwd = 1.5
  )
  abline(
    v = seq(0.5, 2048.5, by = block_length), lty = "dashed",
    lwd = 1.5
  )
  axis(1,
    at = seq(1, length(x_cboe), length.out = 17),
    labels = substr(names(x_cboe)[seq(1, length(x_cboe), length.out = 17)], 12, 16)
  )
  dev.off()
}

{
  png(file = file.path(output_path, "CBOEApplication2.png"), width = 1800, height = 900, res = 210)
  par(mar = c(4, 4, 2, 2))
  plot(cboe_selection[1:2048, 2],
    type = "l", xaxt = "n", xlab = "Horário (hh:mm)",
    ylab = "Valor (dólares)",
    lwd = 1.5
  )
  abline(
    v = seq(0.5, 2048.5, by = block_length), lty = "dashed",
    lwd = 1.5
  )
  axis(1,
    at = seq(1, length(x_cboe), length.out = 17),
    labels = substr(names(x_cboe)[seq(1, length(x_cboe), length.out = 17)], 12, 16)
  )
  axis(2, labels = FALSE)
  dev.off()
}


save.image(file.path(workspace_dir, "11_Changepoint_sim_aplicacoes_part_1.RData"))
set.seed(0)
results_cboe <- changepointdetection_fun(x_cboe, block_length, alpha, jmax, b)

################################################################################

#' Format changepoint results into a consolidated LaTeX table
#'
#' @param results List of detection results from changepointdetection_fun
#' @return A dataframe structured for LaTeX export
latex_fun <- function(results) {
  temp <- results

  for (i in 1:4) {
    if (!is.null(temp[[i]])) {
      temp[[i]] <- cbind(temp[[i]], rep(i, nrow(temp[[i]])))
    }
  }

  # return(temp)

  temp <- arrange(
    distinct(as.data.frame(rbind(
      temp[[1]], temp[[2]],
      temp[[3]], temp[[4]]
    ))),
    V2, V1, V3
  )

  temp <- cbind(temp, NA, NA, NA, NA)

  for (i in 1:nrow(temp)) {
    temp[i, 3 + temp[i, 3]] <- "sim"
  }

  for (i in nrow(temp):2) {
    for (j in (nrow(temp) - 1):1) {
      if (temp[i, 1] == temp[j, 1] & temp[i, 2] == temp[j, 2]) {
        temp[j, 3 + temp[i, 3]] <- "sim"
      }
    }
  }

  for (i in nrow(temp):2) {
    for (j in (i - 1):1) {
      if (temp[i, 1] == temp[j, 1] & temp[i, 2] == temp[j, 2]) {
        print(c(i, j))

        temp <- temp[-i, ]

        break
      }
    }
  }

  return(temp)
}

latex_oceanshear <- latex_fun(results_oceanshear)[, -3]
latex_cboe <- latex_fun(results_cboe)[, -3]

save.image(file.path(workspace_dir, "11_Changepoint_sim_aplicacoes_final.RData"))
