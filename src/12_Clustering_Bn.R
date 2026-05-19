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

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "1_Simulation_functions.R"))
source(file.path(base_path, "src", "2_Bootstrap_methods.R"))
source(file.path(base_path, "src", "7_Quasi_U_statistics_functions.R"))

# --- Testing Mode ---
test_mode <- TRUE
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
bn_matrices_list <- function(data_wv, n, clusters_list) {
  output <- list()

  for (nprime in 2:floor(n / 2)) {
    quasi_u_statistics <- matrix(NA, nrow = nrow(clusters_list[[1]][[nprime - 1]]), ncol = ncol(data_wv))


    for (i in 1:nrow(quasi_u_statistics)) {
      for (j in 1:ncol(quasi_u_statistics)) {
        quasi_u_statistics[i, j] <- b_n(
          data_wv[clusters_list[[1]][[nprime - 1]][i, ], j],
          data_wv[clusters_list[[2]][[nprime - 1]][i, ], j],
          kernel
        )
      }
    }

    output[[nprime - 1]] <- quasi_u_statistics
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

  gen_fun <- get(paste0("model_", model, "_sim"))

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
#' @param x Data matrix
#' @return Matrix with [max_row, col_index, max_value]
max_per_col <- function(x) {
  y <- matrix(0, nrow = ncol(x), ncol = 3)

  for (i in 1:ncol(x)) {
    # Find the index of the maximum element in the current column
    max_index <- which.max(x[, i])

    # Insert values into matrix y
    y[i, 1] <- max_index
    y[i, 2] <- i
    y[i, 3] <- x[max_index, i]
  }

  return(y)
}

#' Apply max_per_col across multiple clusters
#'
#' @param bn_matrices List of Bn matrices
#' @param n Number of elements
#' @return List of results from max_per_col
max_per_col_list <- function(bn_matrices, n) {
  max_bn_scales <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    max_bn_scales[[nprime - 1]] <- max_per_col(bn_matrices[[nprime - 1]])
  }

  return(max_bn_scales)
}

#' Bootstrap Bn statistic across scales
#'
#' @param data_wv_list List containing wavelet variances for two groups
#' @return A vector of pooled bootstrap Bn estimates
bootstrap_bn_multiscale <- function(data_wv_list) {
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
    output <- c(output, b_n(
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
#' @param b Number of bootstrap reps
#' @return List of bootstrap Bn matrices
bootstrap_bn_multiscale_list <- function(data_wv, clusters_list, n, b) {
  bootstrap_bn_versions <- vector(mode = "list", length = floor(n / 2) - 1)

  for (b in 1:b) {
    for (nprime in 2:floor(n / 2)) {
      data_wv_list <- list(
        data_wv[clusters_list[[1]][[nprime - 1]][1, ], ],
        data_wv[clusters_list[[2]][[nprime - 1]][1, ], ]
      )

      bootstrap_bn_versions[[nprime - 1]] <- rbind(
        bootstrap_bn_versions[[nprime - 1]],
        bootstrap_bn_multiscale(data_wv_list)
      )
    }
  }

  return(bootstrap_bn_versions)
}

#' Estimate standard errors from bootstrap samples
#'
#' @param bootstrap_bn_versions List of bootstrap Bn matrices
#' @param n Total elements
#' @return A list of standard error vectors
std_errors_calc <- function(bootstrap_bn_versions, n) {
  std_errors <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    std_errors[[nprime - 1]] <- apply(bootstrap_bn_versions[[nprime - 1]], 2, function(x) {
      sd(x)
    })
  }

  return(std_errors)
}

#' Calculate p-values for observed Bn maxima
#'
#' @param max_bn_scales List of observed maxima
#' @param std_errors List of estimated standard errors
#' @param n Total elements
#' @return List of p-values
p_values_calc <- function(max_bn_scales, std_errors, n) {
  p_values <- vector(mode = "list", length = floor(n / 2) - 1)

  for (nprime in 2:floor(n / 2)) {
    for (j in 1:nrow(max_bn_scales[[1]])) {
      p_values[[nprime - 1]] <- c(p_values[[nprime - 1]], 1 - pnorm(max_bn_scales[[nprime - 1]][j, 3], sd = std_errors[[nprime - 1]][j]))
    }
  }

  return(p_values)
}

#' Identify the best cluster based on minimum p-value
#'
#' @param p_values List of p-values
#' @param t Time series length
#' @param max_bn_scales observed maxima
#' @param clusters_list Possible clusters
#' @return List of best clusters found per scale
best_clusters_fun <- function(p_values, t, max_bn_scales, clusters_list) {
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
      clusters_list[[1]][[min_index[1]]][max_bn_scales[[min_index[[1]]]][min_index[[2]], 1], ],
      clusters_list[[2]][[min_index[1]]][max_bn_scales[[min_index[[1]]]][min_index[[2]], 1], ]
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

start_time <- Sys.time()

n <- 12

n1 <- n / 2

t <- 256

true_cluster <- c(rep(1, n1), rep(2, n - n1))

data <- data_matrix_gen(t, n, n1, "A")

data_wv <- t(apply(data, 1, wv_estimates))

clusters_list <- possible_clusters(n)

bn_matrices <- bn_matrices_list(data_wv, n, clusters_list)

max_bn_scales <- max_per_col_list(bn_matrices, n)

b <- if (test_mode) 2 else 100

bootstrap_bn_versions <- bootstrap_bn_multiscale_list(data_wv, clusters_list, n, b)

std_errors <- std_errors_calc(bootstrap_bn_versions, n)

p_values <- p_values_calc(max_bn_scales, std_errors, n)

best_clusters <- best_clusters_fun(p_values, t, max_bn_scales, clusters_list)

best_clusters <- lapply(best_clusters, function(x) {
  cluster_notation_conv(x, n)
})

unlist(lapply(best_clusters, function(x) {
  adj.rand.index(x, true_cluster)
}))

#

# ACF, COR, EUCL, PER

dist_obj <- diss(data, "ACF")

clust_acf <- pam(dist_obj, k = 2, diss = TRUE)

adj.rand.index(clust_acf$clustering, true_cluster)

end_time <- Sys.time()

end_time - start_time

# 31.59271 secs (16)

################################################################################

#' Run clustering simulation study
#'
#' @param iterations Simulation reps
#' @param b Bootstrap reps
#' @param n Total elements
#' @param n1 Elements in group 1
#' @param t Time series length
#' @param model Model identifier
#' @return List of ARI results for various methods (wv_qu, ACF, COR, EUCL, PER)
main_fun <- function(iterations, b, n, n1, t, model) {
  output_wv_qu <- NULL
  output_acf <- NULL
  output_cor <- NULL
  output_eucl <- NULL
  output_per <- NULL

  true_cluster <- c(rep(1, n1), rep(2, n - n1))

  for (i in 1:iterations) {
    data <- data_matrix_gen(t, n, n1, "A")

    data_wv <- t(apply(data, 1, wv_estimates))

    clusters_list <- possible_clusters(n)

    bn_matrices <- bn_matrices_list(data_wv, n, clusters_list)

    max_bn_scales <- max_per_col_list(bn_matrices, n)

    bootstrap_bn_versions <- bootstrap_bn_multiscale_list(data_wv, clusters_list, n, b)

    std_errors <- std_errors_calc(bootstrap_bn_versions, n)

    p_values <- p_values_calc(max_bn_scales, std_errors, n)

    best_clusters <- best_clusters_fun(p_values, t, max_bn_scales, clusters_list)

    best_clusters <- lapply(best_clusters, function(x) {
      cluster_notation_conv(x, n)
    })

    output_wv_qu <- rbind(
      output_wv_qu,
      unlist(lapply(best_clusters, function(x) {
        adj.rand.index(x, true_cluster)
      }))
    )

    # ACF dissimilarity

    dist_obj_acf <- diss(data, "ACF")

    clust_acf <- pam(dist_obj_acf, k = 2, diss = TRUE)

    output_acf <- c(output_acf, adj.rand.index(clust_acf$clustering, true_cluster))

    # COR dissimilarity

    dist_obj_cor <- diss(data, "COR")

    clust_cor <- pam(dist_obj_cor, k = 2, diss = TRUE)

    output_cor <- c(output_cor, adj.rand.index(clust_cor$clustering, true_cluster))

    # EUCL dissimilarity

    dist_obj_eucl <- diss(data, "EUCL")

    clust_eucl <- pam(dist_obj_eucl, k = 2, diss = TRUE)

    output_eucl <- c(output_eucl, adj.rand.index(clust_eucl$clustering, true_cluster))

    # PER dissimilarity

    dist_obj_per <- diss(data, "PER")

    clust_per <- pam(dist_obj_per, k = 2, diss = TRUE)

    output_per <- c(output_per, adj.rand.index(clust_per$clustering, true_cluster))
  }

  return(list(output_wv_qu, output_acf, output_cor, output_eucl, output_per))
}

n <- 16

n1 <- n / 2

t <- 512

b <- if (test_mode) 2 else 100

iterations <- if (test_mode) 2 else 5

model <- "A"

start_time <- Sys.time()

test <- main_fun(iterations, b, n, n1, t, model)

end_time <- Sys.time()

end_time - start_time

################################################################################

b <- if (test_mode) 2 else 100

iterations <- if (test_mode) 2 else 100

################################################################################

# Time series length = 128, sample size = 12, n1/n2 approx 1

set.seed(72)

start_time <- Sys.time()

n <- 12

n1 <- n / 2

t <- 128

results_128_12_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_12_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_12_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_12_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 128, sample size = 12, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_1.RData"))
set.seed(73)

start_time <- Sys.time()

n <- 12

n1 <- floor((2 / 3) * n)

t <- 128

results_128_12_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_12_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_12_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_12_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 128, sample size = 14, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_2.RData"))
set.seed(74)

start_time <- Sys.time()

n <- 14

n1 <- n / 2

t <- 128

results_128_14_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_14_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_14_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_14_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 128, sample size = 14, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_3.RData"))
set.seed(75)

start_time <- Sys.time()

n <- 14

n1 <- floor((2 / 3) * n)

t <- 128

results_128_14_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_14_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_14_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_14_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 128, sample size = 16, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_4.RData"))
set.seed(76)

start_time <- Sys.time()

n <- 16

n1 <- n / 2

t <- 128

results_128_16_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_16_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_16_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_16_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 128, sample size = 16, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_5.RData"))
set.seed(77)

start_time <- Sys.time()

n <- 16

n1 <- floor((2 / 3) * n)

t <- 128

results_128_16_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_128_16_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_128_16_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_128_16_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

################################################################################

# Time series length = 256, sample size = 12, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_6.RData"))
set.seed(78)

start_time <- Sys.time()

n <- 12

n1 <- n / 2

t <- 256

results_256_12_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_12_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_12_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_12_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 256, sample size = 12, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_7.RData"))
set.seed(79)

start_time <- Sys.time()

n <- 12

n1 <- floor((2 / 3) * n)

t <- 256

results_256_12_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_12_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_12_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_12_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 256, sample size = 14, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_8.RData"))
set.seed(80)

start_time <- Sys.time()

n <- 14

n1 <- n / 2

t <- 256

results_256_14_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_14_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_14_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_14_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 256, sample size = 14, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_9.RData"))
set.seed(81)

start_time <- Sys.time()

n <- 14

n1 <- floor((2 / 3) * n)

t <- 256

results_256_14_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_14_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_14_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_14_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 256, sample size = 16, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_10.RData"))
set.seed(82)

start_time <- Sys.time()

n <- 16

n1 <- n / 2

t <- 256

results_256_16_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_16_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_16_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_16_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 256, sample size = 16, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_11.RData"))
set.seed(83)

start_time <- Sys.time()

n <- 16

n1 <- floor((2 / 3) * n)

t <- 256

results_256_16_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_256_16_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_256_16_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_256_16_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time


################################################################################

# Time series length = 512, sample size = 12, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_12.RData"))
set.seed(84)

start_time <- Sys.time()

n <- 12

n1 <- n / 2

t <- 512

results_512_12_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_12_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_12_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_12_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 512, sample size = 12, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_13.RData"))
set.seed(85)

start_time <- Sys.time()

n <- 12

n1 <- floor((2 / 3) * n)

t <- 512

results_512_12_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_12_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_12_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_12_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 512, sample size = 14, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_14.RData"))
set.seed(86)

start_time <- Sys.time()

n <- 14

n1 <- n / 2

t <- 512

results_512_14_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_14_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_14_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_14_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 512, sample size = 14, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_15.RData"))
set.seed(87)

start_time <- Sys.time()

n <- 14

n1 <- floor((2 / 3) * n)

t <- 512

results_512_14_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_14_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_14_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_14_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 512, sample size = 16, n1/n2 approx 1

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_16.RData"))
set.seed(88)

start_time <- Sys.time()

n <- 16

n1 <- n / 2

t <- 512

results_512_16_1_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_16_1_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_16_1_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_16_1_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time

# Time series length = 512, sample size = 16, n1/n2 approx 2

save.image(file.path(workspace_dir, "12_Clustering_Bn_part_17.RData"))
set.seed(89)

start_time <- Sys.time()

n <- 16

n1 <- floor((2 / 3) * n)

t <- 512

results_512_16_2_a <- main_fun(iterations, b, n, n1, t, "A")

results_512_16_2_b <- main_fun(iterations, b, n, n1, t, "B")

results_512_16_2_c <- main_fun(iterations, b, n, n1, t, "C")

results_512_16_2_d <- main_fun(iterations, b, n, n1, t, "D")

end_time <- Sys.time()

end_time - start_time


################################################################################

quantiles_fun <- function(list_obj, rounding) {
  return(round(cbind(
    matrix(unlist(lapply(list_obj[2:5], quantile)), nrow = 5),
    apply(list_obj[[1]], 2, quantile)
  ), rounding))
}

quantiles_fun(results_128_12_1_a, 2)

################################################################################

latex_out_fun <- function(list_obj12, list_obj14, list_obj16, rounding, t, q, mod) {
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

latex_out_fun(results_128_12_1_a, results_128_14_1_a, results_128_16_1_a, 2, 128, 1, "A")

latex_out_fun(results_512_12_1_a, results_512_14_1_a, results_512_16_1_a, 2, 128, 1, "A")

for (t in c(128, 256, 512)) {
  for (q in c(1, 2)) {
    for (mod in c("A", "B", "C", "D")) {
      name_obj12 <- paste0("results_", t, "_12_", q, "_", mod)
      name_obj14 <- paste0("results_", t, "_14_", q, "_", mod)
      name_obj16 <- paste0("results_", t, "_16_", q, "_", mod)

      latex_out_fun(get(name_obj12), get(name_obj14), get(name_obj16), 2, t, q, mod)

      cat("\n")
    }
  }
}

################################################################################

# Idea for future work: a heuristic to randomly search the very large (in n) search space.
# It could be pure random search with a time constraint or something based on swapping
# pairs of elements from distinct groups in order to decrease the p-value


save.image(file.path(workspace_dir, "12_Clustering_Bn_final.RData"))
