# =============================================================================
# 0_plots_dados.R
# =============================================================================
# Purpose  : General data visualization for thesis datasets and methods.
# Chapter  : Multiple
# Inputs   : Ocean Shear, SPX, MJO, Arrowhead, INMET (Data files not included).
# Outputs  : Data overview plots (PNGs).
# Depends  : 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running

source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

library(grDevices)

# Set and create output directory for plots
OUTPUT_PATH <- file.path(BASE_PATH, "Plots/Plots_0")
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)


################################################################################

# Cisalhamento Oceânico

# Note: Data file is not included in the repository. See thesis for sources.
ocean_shear <- unname(unlist(as.vector(read.table(file.path(BASE_PATH, "Dados", "Ocean Shear", "Ocean Shear.txt")))))


{
     png(file = file.path(OUTPUT_PATH, "OceanShear.png"), width = 1800, height = 900, res = 210)

     par(mfrow = c(1, 1), mar = c(4.1, 4, 1, 1))
     plot(ocean_shear,
          type = "l", xaxt = "n", xlab = "Profundidade (metros)",
          ylab = "Cisalhamento", yaxt = "n",
          lwd = 1.5
     )
     axis(1,
          at = seq(1, 6875, length.out = 10),
          labels = round(seq(350, 1037.4, length.out = 10))
     )
     axis(2)

     dev.off()
}

oceanshear_diff <- diff(ocean_shear)

{
     png(file = file.path(OUTPUT_PATH, "OceanShear2.png"), width = 1800, height = 900, res = 210)
     par(mfrow = c(1, 1), mar = c(4.1, 5, 1, 1))
     plot(oceanshear_diff,
          type = "l", xaxt = "n", xlab = "Profundidade (metros)",
          ylab = "Cisalhamento \n (diferenciado)", yaxt = "n",
          lwd = 1.5
     )
     axis(1,
          at = seq(2, 6875, length.out = 10),
          labels = round(seq(350.1, 1037.4, length.out = 10))
     )
     axis(2)

     dev.off()
}


N <- length(oceanshear_diff)

n_levels <- floor(log2(1 + (N - 1) / (8 - 1)))

oceanshear_diff.modwt <- modwt(oceanshear_diff, n.levels = n_levels)[1:n_levels]

for (j in 1:n_levels) {
     L_j <- (2^j - 1) * (8 - 1) + 1

     oceanshear_diff.modwt[[j]] <- c(rep(NA, L_j), oceanshear_diff.modwt[[j]][L_j:length(oceanshear_diff.modwt[[j]])])
}

{
     png(file = file.path(OUTPUT_PATH, "OceanShear3.png"), width = 1600, height = 1600, res = 280)

     par(mar = c(3, 3.5, 0.1, 1))
     par(mfrow = c(6, 1))
     par(oma = c(0, 0, 0, 0))

     plot(oceanshear_diff.modwt[[1]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1, labels = FALSE, at = seq(2, 6875, length.out = 10))
     title(ylab = expression(W[list(1, t)]), mgp = c(2, 1, 0))

     plot(oceanshear_diff.modwt[[2]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1, labels = FALSE, at = seq(2, 6875, length.out = 10))
     title(ylab = expression(W[list(2, t)]), mgp = c(2, 1, 0))

     plot(oceanshear_diff.modwt[[3]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1, labels = FALSE, at = seq(2, 6875, length.out = 10))
     title(ylab = expression(W[list(3, t)]), mgp = c(2, 1, 0))

     plot(oceanshear_diff.modwt[[4]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1, labels = FALSE, at = seq(2, 6875, length.out = 10))
     title(ylab = expression(W[list(4, t)]), mgp = c(2, 1, 0))

     plot(oceanshear_diff.modwt[[5]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1, labels = FALSE, at = seq(2, 6875, length.out = 10))
     title(ylab = expression(W[list(5, t)]), mgp = c(2, 1, 0))

     plot(oceanshear_diff.modwt[[6]],
          type = "l", xaxt = "n", yaxt = "n", ylab = "",
          lwd = 1.5
     )
     axis(1,
          at = seq(2, 6875, length.out = 10),
          labels = round(seq(350.1, 1037.4, length.out = 10))
     )
     title(
          ylab = expression(W[list(6, t)]), mgp = c(2, 1, 0),
          xlab = "Profundidade (metros)"
     )

     dev.off()
}

################################################################################

# Fixação de Membrana

################################################################################

# Índice da bolsa americana

#' Calculate relative returns from price data
#'
#' @param data Dataframe with prices (first column dates/times)
#' @return Vector of relative returns
returns_fun <- function(data) {
     data_temp <- unname(unlist(as.vector(data[, 2:ncol(data)])))

     output <- NULL

     for (i in 2:length(data_temp)) {
          output <- c(output, (data_temp[i] - data_temp[i - 1]) / data_temp[i - 1])
     }

     return(output)
}

#' Subset data by time intervals
#'
#' @param interval_type Either "minutes", "hours", or "days"
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
     }

     output1 <- seq(1, 24600, by = adj_interval)

     return(list(output1, c(1, day_start:day_stop)))
}

# Note: Data file is not included in the repository. See thesis for sources.
CBOE <- read.csv(file.path(BASE_PATH, "Dados", "SPX_second", "SPX.csv"), check.names = FALSE)

# sampling
indexes <- subset_data("minutes", 15, 2, 44)

CBOE_selection <- CBOE[indexes[[1]], indexes[[2]]]

CBOE_selection_temp <- unname(unlist(as.vector(CBOE_selection[, 2:ncol(CBOE_selection)])))

labels_plot <- substr(colnames(CBOE_selection)[1 + seq(1, 43, by = 7)], 6, 10)

labels_plot <- sapply(labels_plot, function(x) {
     paste0(substr(x, 4, 5), ".", substr(x, 1, 2))
})

{
     png(file = file.path(OUTPUT_PATH, "CBOE.png"), width = 1800, height = 900, res = 210)

     par(mar = c(4, 4, 1, 1), mfrow = c(1, 1))

     plot(CBOE_selection_temp,
          type = "l", xlab = "Data (dia.mês)",
          ylab = "", xaxt = "n",
          lwd = 1.5
     )
     axis(2, labels = FALSE)
     axis(1,
          at = seq(1, 43 * 28, by = 7 * 28),
          labels = labels_plot
     )
     title(ylab = "Valor (dólares)", mgp = c(3, 1, 0))
     dev.off()
}

CBOE_selection_returns <- returns_fun(CBOE_selection)

{
     png(file = file.path(OUTPUT_PATH, "CBOE2.png"), width = 1800, height = 900, res = 210)

     par(mar = c(4, 4, 1, 1))

     plot(CBOE_selection_returns,
          type = "l", xlab = "Data (dia.mês)",
          xaxt = "n", ylab = "Retorno",
          lwd = 1.5
     )
     axis(1,
          at = seq(1, 43 * 28, by = 7 * 28),
          labels = labels_plot
     )

     dev.off()
}

################################################################################

# Oscilação de Madden-Julian

# Note: Data file is not included in the repository. See thesis for sources.
# Added na.strings="*****" to handle missing values and prevent coercion warnings
MJO_df <- read.table(file.path(BASE_PATH, "Dados", "MJO", "MJO.txt"),
     na.strings = c("NA", "*****")
)
MJO_labels <- substr(rownames(MJO_df)[2:nrow(MJO_df)], 1, 4)

# Data is now numeric automatically thanks to na.strings
MJO <- MJO_df$INDEX_1[2:nrow(MJO_df)]


MJO_paper <- MJO[1:2354]
MJO_labels_paper <- MJO_labels[1:2354]

years_plot <- seq(1978, 2010, by = 2)
indexes_years_plot <- NULL

for (year in years_plot) {
     indexes_years_plot <- c(indexes_years_plot, which(MJO_labels_paper == year)[1])
}

{
     png(file = file.path(OUTPUT_PATH, "MJO.png"), width = 1800, height = 900, res = 210)

     par(mar = c(4.1, 5, 1, 1))
     par(mfrow = c(1, 1))
     plot(MJO_paper,
          type = "l", xaxt = "n", ylim = c(-4, 4), xlab = "Ano",
          ylab = "Potencial de velocidade \n (transformado)",
          lwd = 1.5
     )
     axis(1, at = indexes_years_plot, labels = years_plot)
     dev.off()
}


################################################################################

# Pontas de flechas

# Note: Data file is not included in the repository. See thesis for sources.
arrowhead <- unname(unlist(read.table(file.path(BASE_PATH, "Dados", "Arrowhead", "ArrowHead.txt"))))

{
     indexes_min <- c(1, 64, 126, 189, 251, 314, 377, 439, 503, 565, 629, 691, 753, 816, 878, 942, 1005, 1068, 1130, 1192, 1255, 1319, 1381, 1443, 1506)

     png(
          file = file.path(OUTPUT_PATH, "arrowhead.png"),
          width = 1200, height = 600, res = 140
     )

     par(mar = c(4, 5, 1, 1))
     par(mfrow = c(2, 1))
     plot(arrowhead[1:753],
          type = "l", xlab = "Índice (Avonlea)", ylab = "Distância",
          lwd = 1.5, ylim = c(-2.5, 3), yaxt = "n", xaxt = "n"
     )
     abline(v = indexes_min, col = "black", lty = 2)
     axis(1, at = indexes_min[c(1, 3, 5, 7, 9, 11, 13)], cex.axis = 1.5, cex.lab = 1.5, )
     axis(2, at = c(-2, 0, 2), cex.axis = 1.5)

     plot(arrowhead[754:1506],
          type = "l", xlab = "Índice (Clovis)", ylab = "Distância",
          lwd = 1.5, ylim = c(-2.5, 3), yaxt = "n", xaxt = "n"
     )
     abline(v = indexes_min, col = "black", lty = 2)
     axis(1, at = indexes_min[c(13, 15, 17, 19, 21, 23, 25)] - 753, labels = indexes_min[c(13, 15, 17, 19, 21, 23, 25)], cex.axis = 1.5, cex.lab = 1.5)
     axis(2, at = c(-2, 0, 2), cex.axis = 1.5)

     dev.off()
}

################################################################################

# Quebra de rigidez dielétrica

################################################################################

# Temperaturas em cidades Brasileiras

# Note: Data file is not included in the repository. See thesis for sources.
INMET <- read.csv(file.path(BASE_PATH, "Dados", "Selecao", "dados_INMET_processados.csv"))

INMET_selec <- INMET[seq(1, 8760, by = 12), ]

first_days <- NULL

for (i in unique(substr(INMET_selec[, 1], 1, 2))) {
     first_days <- c(first_days, which(substr(INMET_selec[, 1], 1, 2) == i)[1])
}

{
     png(
          file = file.path(OUTPUT_PATH, "Temperaturas.png"),
          width = 1800, height = 1200, res = 210
     )

     par(mar = c(4.1, 4.1, 1, 1))
     par(mfrow = c(1, 1))
     plot(INMET_selec[, 3],
          type = "l", ylim = c(9, 35),
          xlab = "Mês", ylab = "Temperatura (°C)",
          xaxt = "n",
          lwd = 1.5
     )
     lines(INMET_selec[, 10], col = "grey")
     # axis(1,at = seq(0,609, length.out = 6), labels = c("Jan", "Mar", "Mai", "Jul", "Set", "Nov"))
     axis(1, at = c(first_days, length(INMET_selec[, 3] + 1)), labels = c(
          "Jan", "Fev", "Mar",
          "Abr", "Mai", "Jun",
          "Jul", "Ago", "Set",
          "Out", "Nov", "Dez", "Jan"
     ))
     legend("bottomleft",
          legend = c("Altamira", "Maringá"),
          col = c("black", "grey"),
          lty = c("solid", "solid"),
          lwd = 1.5
     )

     dev.off()
}

INMET_selec_diff <- apply(INMET_selec[-(1:2)], 2, diff)

{
     png(
          file = file.path(OUTPUT_PATH, "Temperaturas2.png"),
          width = 1800, height = 1200, res = 210
     )

     par(mfrow = c(2, 1))
     par(mar = c(3, 5, 1, 1))
     plot(INMET_selec_diff[, 1],
          type = "l", ylim = c(-14, 14),
          xlab = "", ylab = "Temperatura diferenciada \n (Altamira)",
          xaxt = "n",
          lwd = 1.5
     )
     axis(1, at = c(first_days, length(INMET_selec[, 3] + 1)), labels = FALSE)

     plot(INMET_selec_diff[, 8],
          type = "l", ylim = c(-14, 14),
          xlab = "", ylab = "Temperatura diferenciada \n (Maringá)", col = "grey",
          xaxt = "n",
          lwd = 1.5
     )
     axis(1, at = c(first_days, length(INMET_selec[, 3] + 1)), labels = c(
          "Jan", "Fev", "Mar",
          "Abr", "Mai", "Jun",
          "Jul", "Ago", "Set",
          "Out", "Nov", "Dez", "Jan"
     ))
     title(xlab = "Mês", mgp = c(2, 1, 0))
     dev.off()
}

################################################################################

# Exemplo Bootstrap

set.seed(0)

boot_example <- 0:127

boot_example_extended <- c(boot_example, boot_example)

block_lengths <- rgeom(6, 1 / 28)

block_start_points <- sample(128, 6, replace = TRUE)

boot_example_bootstrap <- c(
     boot_example_extended[block_start_points[1]:(block_start_points[1] + block_lengths[1] - 1)],
     boot_example_extended[block_start_points[2]:(block_start_points[2] + block_lengths[2] - 1)],
     boot_example_extended[block_start_points[3]:(block_start_points[3] + block_lengths[3] - 1)],
     boot_example_extended[block_start_points[4]:(block_start_points[4] + block_lengths[4] - 1)],
     boot_example_extended[block_start_points[5]:(block_start_points[5] + block_lengths[5] - 1)],
     boot_example_extended[block_start_points[6]:(block_start_points[6] + block_lengths[6] - 1)]
)

{
     png(file = file.path(OUTPUT_PATH, "boot_example_1.png"), width = 1800, height = 900, res = 210)

     par(mar = c(3.1, 3, 1, 1))
     par(mfrow = c(1, 1))
     plot(boot_example, type = "p", pch = 1, xlab = "", ylab = "", lwd = 1.5)
     lines(boot_example, lwd = 1.5)
     dev.off()
}

{
     png(file = file.path(OUTPUT_PATH, "boot_example_2.png"), width = 1800, height = 900, res = 210)

     par(mar = c(3.1, 3, 1, 1))
     par(mfrow = c(1, 1))
     plot(boot_example_bootstrap, type = "p", pch = 1, xlab = "", ylab = "", lwd = 1.5)
     lines(boot_example_bootstrap, lwd = 1.5)
     abline(v = c(5.5, 12.5, 22.5, 97.5, 116.5), lty = 2)
     dev.off()
}
