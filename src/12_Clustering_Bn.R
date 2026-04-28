# =============================================================================
# 12_Clustering_Bn.R
# =============================================================================
# Purpose  : Simulation study for time series clustering based on Bn statistics.
# Chapter  : Chapter 4, Appendix C.2
# Inputs   : 1_Simulation_functions.R, 2_Bootstrap_methods.R, 7_Quasi_U_statistics_functions.R
# Outputs  : ARI comparisons against ACF, PER, COR, and EUCL Dissimilarities.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R, 7_Quasi_U_statistics_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
WORKSPACE_DIR <- file.path(BASE_PATH, "src", "WorkspaceData")
if (!dir.exists(WORKSPACE_DIR)) dir.create(WORKSPACE_DIR, recursive = TRUE)

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))
source(file.path(BASE_PATH, "src", "7_Quasi_U_statistics_functions.R"))

# --- Testing Mode ---
TEST_MODE <- TRUE
# --------------------

# Define kernel for Bn
kernel <- function(x, y) (x - y)^2

library(fossil)

library(TSclust)

library(cluster)

################################################################################

#' Calculate all possible two-group clusters
#'
#' @param n Number of elements
#' @return A list with two lists: clusters_1 and clusters_2
possible_clusters <- function(n) {
  output1 <- list()
  output2 <- list()

  for (nprime in 2:floor(n / 2)) {
    clusters_1 <- t(combn(n, nprime))

    if (nprime == n / 2) {
      clusters_1 <- clusters_1[1:(nrow(clusters_1) / 2), ]
    }

    output1[[nprime - 1]] <- clusters_1

    clusters_2 <- t(apply(clusters_1, 1, function(x) {
      setdiff(1:n, x)
    }))

    output2[[nprime - 1]] <- clusters_2
  }

  return(list(output1, output2))
}

#' Calculate Bn matrices for all clusters across scales
#'
#' @param data_wv Matrix of wavelet variance estimates
#' @param n Number of elements
#' @param clusters_list Output from possible_clusters
#' @return A list containing Bn matrices for each size division
Bn_matrices_list <- function(data_wv, n, clusters_list) {
  output <- list()

  for (nprime in 2:floor(n / 2)) {
    quasi_U_statistics <- matrix(NA, nrow = nrow(clusters_list[[1]][[nprime - 1]]), ncol = ncol(data_wv))


    for (i in 1:nrow(quasi_U_statistics)) {
      for (j in 1:ncol(quasi_U_statistics)) {
        quasi_U_statistics[i, j] <- B_n(
          data_wv[clusters_list[[1]][[nprime - 1]][i, ], j],
          data_wv[clusters_list[[2]][[nprime - 1]][i, ], j],
          kernel
        )
      }
    }

    output[[nprime - 1]] <- quasi_U_statistics
  }

  return(output)
}

#' Generate simulated data matrix for clustering
#'
#' @param t Time series length
#' @param n Total number of time series
#' @param n1 Number of series in first group
#' @param model Model identifier ("A" through "D")
#' @return A matrix with simulated time series
data_matrix_gen <- function(t, n, n1, model) {
  data <- matrix(NA, nrow = n, ncol = t)

  gen_fun <- get(paste0("Model_", model, "_sim"))

  for (i in 1:n) {
    if (i <= n1) {
      data[i, ] <- gen_fun(t)
    } else {
      data[i, ] <- 1.25 * gen_fun(t)
    }
  }

  return(data)
}

#' Find maximum value and index per column
#'
#' @param X Data matrix
#' @return Matrix with [max_row, col_index, max_value]
max_per_col <- function(X) {
  Y <- matrix(0, nrow = ncol(X), ncol = 3)

  for (i in 1:ncol(X)) {
    # Find the index of the maximum element in the current column
    max_index <- which.max(X[, i])

    # Insert values into matrix Y
    Y[i, 1] <- max_index
    Y[i, 2] <- i
    Y[i, 3] <- X[max_index, i]
  }

  return(Y)
}

#' Apply max_per_col across multiple clusters
#'
#' @param Bn_matrices List of Bn matrices
#' @param n Number of elements
#' @return List of results from max_per_col
max_per_col_list <- function(Bn_matrices, n) {
  max_Bn_scales <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    max_Bn_scales[[nprime - 1]] <- max_per_col(Bn_matrices[[nprime - 1]])
  }

  return(max_Bn_scales)
}

#' Bootstrap Bn statistic across scales
#'
#' @param data_wv_list List containing wavelet variances for two groups
#' @return A vector of pooled bootstrap Bn estimates
Bootstrap_Bn_multiscale <- function(data_wv_list) {
  data1 <- data_wv_list[[1]]
  data2 <- data_wv_list[[2]]

  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n <- n1 + n2

  data <- rbind(data1, data2)

  indexes <- sample(n, replace = TRUE)

  data <- data[indexes, ]

  data1_temp <- data[1:n1, ]
  data2_temp <- data[(n1 + 1):n, ]

  # return(list(data1_temp, data2_temp))

  output <- NULL

  for (j in 1:ncol(data1)) {
    output <- c(output, B_n(
      data1_temp[, j],
      data2_temp[, j],
      kernel
    ))
  }

  return(output)
}

#' Generate bootstrap Bn versions for all cluster size divisions
#'
#' @param data_wv Matrix of wavelet variance estimates
#' @param clusters_list List of possible clusters
#' @param n Total elements
#' @param B Number of bootstrap reps
#' @return List of bootstrap Bn matrices
Bootstrap_Bn_multiscale_list <- function(data_wv, clusters_list, n, B) {
  Bootstrap_Bn_versions <- vector(mode = "list", length = floor(n / 2) - 1)

  for (b in 1:B) {
    for (nprime in 2:floor(n / 2)) {
      data_wv_list <- list(
        data_wv[clusters_list[[1]][[nprime - 1]][1, ], ],
        data_wv[clusters_list[[2]][[nprime - 1]][1, ], ]
      )

      Bootstrap_Bn_versions[[nprime - 1]] <- rbind(
        Bootstrap_Bn_versions[[nprime - 1]],
        Bootstrap_Bn_multiscale(data_wv_list)
      )
    }
  }

  return(Bootstrap_Bn_versions)
}

#' Estimate standard errors from bootstrap samples
#'
#' @param Bootstrap_Bn_versions List of bootstrap Bn matrices
#' @param n Total elements
#' @return A list of standard error vectors
std.errors_calc <- function(Bootstrap_Bn_versions, n) {
  std.errors <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    std.errors[[nprime - 1]] <- apply(Bootstrap_Bn_versions[[nprime - 1]], 2, function(x) {
      sd(x)
    })
  }

  return(std.errors)
}

#' Calculate p-values for observed Bn maxima
#'
#' @param max_Bn_scales List of observed maxima
#' @param std.errors List of estimated standard errors
#' @param n Total elements
#' @return List of p-values
p_values_calc <- function(max_Bn_scales, std.errors, n) {
  p_values <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    for (j in 1:nrow(max_Bn_scales[[1]])) {
      p_values[[nprime - 1]] <- c(p_values[[nprime - 1]], 1 - pnorm(max_Bn_scales[[nprime - 1]][j, 3], sd = std.errors[[nprime - 1]][j]))
    }
  }

  return(p_values)
}

#' Identify the best cluster based on minimum p-value
#'
#' @param p_values List of p-values
#' @param t Time series length
#' @param max_Bn_scales observed maxima
#' @param clusters_list Possible clusters
#' @return List of best clusters found per scale
best_clusters_fun <- function(p_values, t, max_Bn_scales, clusters_list) {
  min_value <- Inf
  min_index <- NULL

  max_scale <- log2(t) - 3

  best_cluster <- vector(mode = "list", length = max_scale)

  for (scal in 1:max_scale) {
    for (i in 1:length(p_values)) {
      current_array <- p_values[[i]][1:scal]

      current_min <- min(current_array)
      current_min_index <- which(current_array == current_min, arr.ind = TRUE)

      if (current_min < min_value) {
        min_value <- current_min
        min_index <- c(i, current_min_index)
      }
    }

    best_cluster[[scal]] <- list(
      clusters_list[[1]][[min_index[1]]][max_Bn_scales[[min_index[[1]]]][min_index[[2]], 1], ],
      clusters_list[[2]][[min_index[1]]][max_Bn_scales[[min_index[[1]]]][min_index[[2]], 1], ]
    )
  }


  return(best_cluster)
}

#' Convert cluster list to membership notation
#'
#' @param cluster_list List of element indexes for each cluster
#' @param n Total elements
#' @return Vector of membership labels
cluster_notation_conv <- function(cluster_list, n) {
  output <- rep(0, n)

  for (i in 1:length(cluster_list)) {
    output[cluster_list[[i]]] <- i
  }

  return(output)
}

################################################################################

{
  start.time <- Sys.time()

  n <- 12

  n1 <- n / 2

  t <- 256

  true_cluster <- c(rep(1, n1), rep(2, n - n1))

  data <- data_matrix_gen(t, n, n1, "A")

  data_wv <- t(apply(data, 1, wv_estimates))

  clusters_list <- possible_clusters(n)

  Bn_matrices <- Bn_matrices_list(data_wv, n, clusters_list)

  max_Bn_scales <- max_per_col_list(Bn_matrices, n)

  B <- if (TEST_MODE) 2 else 100

  Bootstrap_Bn_versions <- Bootstrap_Bn_multiscale_list(data_wv, clusters_list, n, B)

  std.errors <- std.errors_calc(Bootstrap_Bn_versions, n)

  p_values <- p_values_calc(max_Bn_scales, std.errors, n)

  best_clusters <- best_clusters_fun(p_values, t, max_Bn_scales, clusters_list)

  best_clusters <- lapply(best_clusters, function(x) {
    cluster_notation_conv(x, n)
  })

  unlist(lapply(best_clusters, function(x) {
    adj.rand.index(x, true_cluster)
  }))

  #

  # ACF, COR, EUCL, PER

  dist_obj <- diss(data, "ACF")

  clust_ACF <- pam(dist_obj, k = 2, diss = TRUE)

  adj.rand.index(clust_ACF$clustering, true_cluster)

  end.time <- Sys.time()

  end.time - start.time
}

# 31.59271 secs (16)

################################################################################

#' Run clustering simulation study
#'
#' @param iterations Simulation reps
#' @param B Bootstrap reps
#' @param n Total elements
#' @param n1 Elements in group 1
#' @param t Time series length
#' @param model Model identifier
#' @return List of ARI results for various methods (wv_qU, ACF, COR, EUCL, PER)
main_fun <- function(iterations, B, n, n1, t, model) {
  output_wv_qU <- NULL
  output_ACF <- NULL
  output_COR <- NULL
  output_EUCL <- NULL
  output_PER <- NULL

  true_cluster <- c(rep(1, n1), rep(2, n - n1))

  for (i in 1:iterations) {
    data <- data_matrix_gen(t, n, n1, "A")

    data_wv <- t(apply(data, 1, wv_estimates))

    clusters_list <- possible_clusters(n)

    Bn_matrices <- Bn_matrices_list(data_wv, n, clusters_list)

    max_Bn_scales <- max_per_col_list(Bn_matrices, n)

    Bootstrap_Bn_versions <- Bootstrap_Bn_multiscale_list(data_wv, clusters_list, n, B)

    std.errors <- std.errors_calc(Bootstrap_Bn_versions, n)

    p_values <- p_values_calc(max_Bn_scales, std.errors, n)

    best_clusters <- best_clusters_fun(p_values, t, max_Bn_scales, clusters_list)

    best_clusters <- lapply(best_clusters, function(x) {
      cluster_notation_conv(x, n)
    })

    output_wv_qU <- rbind(
      output_wv_qU,
      unlist(lapply(best_clusters, function(x) {
        adj.rand.index(x, true_cluster)
      }))
    )

    # ACF dissimilarity

    dist_obj_ACF <- diss(data, "ACF")

    clust_ACF <- pam(dist_obj_ACF, k = 2, diss = TRUE)

    output_ACF <- c(output_ACF, adj.rand.index(clust_ACF$clustering, true_cluster))

    # COR dissimilarity

    dist_obj_COR <- diss(data, "COR")

    clust_COR <- pam(dist_obj_COR, k = 2, diss = TRUE)

    output_COR <- c(output_COR, adj.rand.index(clust_COR$clustering, true_cluster))

    # EUCL dissimilarity

    dist_obj_EUCL <- diss(data, "EUCL")

    clust_EUCL <- pam(dist_obj_EUCL, k = 2, diss = TRUE)

    output_EUCL <- c(output_EUCL, adj.rand.index(clust_EUCL$clustering, true_cluster))

    # PER dissimilarity

    dist_obj_PER <- diss(data, "PER")

    clust_PER <- pam(dist_obj_PER, k = 2, diss = TRUE)

    output_PER <- c(output_PER, adj.rand.index(clust_PER$clustering, true_cluster))
  }

  return(list(output_wv_qU, output_ACF, output_COR, output_EUCL, output_PER))
}

n <- 16

n1 <- n / 2

t <- 512

B <- if (TEST_MODE) 2 else 100

iterations <- if (TEST_MODE) 2 else 5

model <- "A"

{
  start.time <- Sys.time()

  test <- main_fun(iterations, B, n, n1, t, model)

  end.time <- Sys.time()

  end.time - start.time
}

################################################################################

B <- if (TEST_MODE) 2 else 100

iterations <- if (TEST_MODE) 2 else 100

################################################################################

# Time series length = 128, sample size = 12, n1/n2 approx 1

set.seed(72)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- n / 2

  t <- 128

  results_128_12_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_12_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_12_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_12_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 128, sample size = 12, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_1.RData"))
set.seed(73)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- floor((2 / 3) * n)

  t <- 128

  results_128_12_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_12_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_12_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_12_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 128, sample size = 14, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_2.RData"))
set.seed(74)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- n / 2

  t <- 128

  results_128_14_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_14_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_14_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_14_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 128, sample size = 14, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_3.RData"))
set.seed(75)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- floor((2 / 3) * n)

  t <- 128

  results_128_14_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_14_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_14_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_14_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 128, sample size = 16, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_4.RData"))
set.seed(76)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- n / 2

  t <- 128

  results_128_16_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_16_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_16_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_16_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 128, sample size = 16, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_5.RData"))
set.seed(77)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- floor((2 / 3) * n)

  t <- 128

  results_128_16_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_128_16_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_128_16_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_128_16_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

################################################################################

# Time series length = 256, sample size = 12, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_6.RData"))
set.seed(78)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- n / 2

  t <- 256

  results_256_12_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_12_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_12_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_12_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 256, sample size = 12, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_7.RData"))
set.seed(79)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- floor((2 / 3) * n)

  t <- 256

  results_256_12_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_12_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_12_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_12_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 256, sample size = 14, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_8.RData"))
set.seed(80)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- n / 2

  t <- 256

  results_256_14_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_14_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_14_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_14_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 256, sample size = 14, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_9.RData"))
set.seed(81)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- floor((2 / 3) * n)

  t <- 256

  results_256_14_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_14_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_14_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_14_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 256, sample size = 16, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_10.RData"))
set.seed(82)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- n / 2

  t <- 256

  results_256_16_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_16_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_16_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_16_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 256, sample size = 16, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_11.RData"))
set.seed(83)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- floor((2 / 3) * n)

  t <- 256

  results_256_16_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_256_16_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_256_16_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_256_16_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}


################################################################################

# Time series length = 512, sample size = 12, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_12.RData"))
set.seed(84)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- n / 2

  t <- 512

  results_512_12_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_12_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_12_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_12_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 512, sample size = 12, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_13.RData"))
set.seed(85)

{
  start.time <- Sys.time()

  n <- 12

  n1 <- floor((2 / 3) * n)

  t <- 512

  results_512_12_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_12_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_12_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_12_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 512, sample size = 14, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_14.RData"))
set.seed(86)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- n / 2

  t <- 512

  results_512_14_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_14_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_14_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_14_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 512, sample size = 14, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_15.RData"))
set.seed(87)

{
  start.time <- Sys.time()

  n <- 14

  n1 <- floor((2 / 3) * n)

  t <- 512

  results_512_14_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_14_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_14_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_14_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 512, sample size = 16, n1/n2 approx 1

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_16.RData"))
set.seed(88)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- n / 2

  t <- 512

  results_512_16_1_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_16_1_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_16_1_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_16_1_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}

# Time series length = 512, sample size = 16, n1/n2 approx 2

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_part_17.RData"))
set.seed(89)

{
  start.time <- Sys.time()

  n <- 16

  n1 <- floor((2 / 3) * n)

  t <- 512

  results_512_16_2_A <- main_fun(iterations, B, n, n1, t, "A")

  results_512_16_2_B <- main_fun(iterations, B, n, n1, t, "B")

  results_512_16_2_C <- main_fun(iterations, B, n, n1, t, "C")

  results_512_16_2_D <- main_fun(iterations, B, n, n1, t, "D")

  end.time <- Sys.time()

  end.time - start.time
}


################################################################################

quantiles_fun <- function(list_obj, rounding) {
  return(round(cbind(
    matrix(unlist(lapply(list_obj[2:5], quantile)), nrow = 5),
    apply(list_obj[[1]], 2, quantile)
  ), rounding))
}

quantiles_fun(results_128_12_1_A, 2)

################################################################################

Latex_out_fun <- function(list_obj12, list_obj14, list_obj16, rounding, t, q, mod) {
  levels <- ncol(list_obj12[[1]])

  cat("\\begin{table}[h!] \n")

  cat("\\centering \n")

  if (q == 1) {
    temp <- "= "
  } else {
    temp <- "\\approx "
  }

  cat(paste0(
    "\\caption{Quartis empíricos dos índices ajustados de Rand para $t = ", t, "$, ",
    "$n_1/n_2 ", temp, q, "$ e modelo (", tolower(mod), ")", ".} \n"
  ))

  # cat(paste0("\\caption{Empirical quantiles of the adjusted Rand indexes for $t = ",  t  ,"$, ",
  # "$n = ", n, "$, " ,"$n_1/n_2 " , temp, q, "$ and model (", tolower(mod), ")", ".} \n"))

  cat(paste0("\\label{tab_clust_", t, "_", q, "_", tolower(mod), "} \n"))

  cat(paste0(
    "\\begin{tabular}{cc|",
    strrep("c", levels + 4),
    "} \n"
  ))

  # Header

  header <- "Quart. & n & ACF & COR & EUCL & PER & \\#j = 1"

  for (i in 2:(levels)) {
    header <- paste0(header, " & $", i, "$")
  }

  header <- paste0(header, " \\\\ \\hline \n")

  cat(header)

  for (i in c(2, 4)) {
    cat(paste0(" & ", 12, " & ", paste0(numformat(quantiles_fun(list_obj12, 2)[i, ]), collapse = " & "), " \\\\ \n"))

    cat(paste0(25 * (i - 1), "\\% & ", 14, " & ", paste0(numformat(quantiles_fun(list_obj14, 2)[i, ]), collapse = " & "), " \\\\ \n"))

    cat(paste0(" & ", 16, " & ", paste0(numformat(quantiles_fun(list_obj16, 2)[i, ]), collapse = " & "), " \\\\ \\hline \n"))
  }

  cat("\\end{tabular} \n")

  cat("\\end{table} \n")
}

Latex_out_fun(results_128_12_1_A, results_128_14_1_A, results_128_16_1_A, 2, 128, 1, "A")

Latex_out_fun(results_512_12_1_A, results_512_14_1_A, results_512_16_1_A, 2, 128, 1, "A")

for (t in c(128, 256, 512)) {
  for (q in c(1, 2)) {
    for (mod in c("A", "B", "C", "D")) {
      name_obj12 <- paste0("results_", t, "_12_", q, "_", mod)
      name_obj14 <- paste0("results_", t, "_14_", q, "_", mod)
      name_obj16 <- paste0("results_", t, "_16_", q, "_", mod)

      Latex_out_fun(get(name_obj12), get(name_obj14), get(name_obj16), 2, t, q, mod)

      cat("\n")
    }
  }
}

################################################################################

# Idea for future work: a heuristic to randomly search the very large (in n) search space.
# It could be pure random search with a time constraint or something based on swapping
# pairs of elements from distinct groups in order to decrease the p-value


save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_final.RData"))
