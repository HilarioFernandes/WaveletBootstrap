# =============================================================================
# 12_Clustering_Bn_aplicacoes.R
# =============================================================================
# Purpose  : Real-data clustering applications (Arrowhead, INMET).
# Chapter  : Chapter 4
# Inputs   : ArrowHead.txt, dados_INMET_processados.csv (Data files not included).
# Outputs  : Clustering labels and ARI results for various methods.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R, 7_Quasi_U_statistics_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running
WORKSPACE_DIR <- file.path(BASE_PATH, "src", "WorkspaceData")
if(!dir.exists(WORKSPACE_DIR)) dir.create(WORKSPACE_DIR, recursive=TRUE)

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

B <- if(TEST_MODE) 5 else 100

################################################################################

#' Calculate all possible two-group clusters
#' 
#' @param n Number of elements
#' @return A list with two lists: clusters_1 and clusters_2
possible_clusters <- function(n){
  
  output1 <- list()
  output2 <- list()
  
  for(nprime in 2:floor(n/2)){
    
    clusters_1 <- t(combn(n,nprime))
    
    if(nprime == n/2){
      
      clusters_1 <- clusters_1[1:(nrow(clusters_1)/2),]
      
    }
    
    output1[[nprime-1]] <- clusters_1
    
    clusters_2 <- t(apply(clusters_1, 1, function(x){setdiff(1:n, x)}))
    
    output2[[nprime-1]] <- clusters_2
    
  }
  
  return(list(output1, output2))
  
}

#' Calculate Bn matrices for all clusters across scales
#' 
#' @param data_wv Matrix of wavelet variance estimates
#' @param n Number of elements
#' @param clusters_list Output from possible_clusters
#' @return A list containing Bn matrices for each size division
Bn_matrices_list <- function(data_wv,n,clusters_list){
  
  output <- list()
  
  for(nprime in 2:floor(n/2)){
    
    quasi_U_statistics <- matrix(NA, nrow = nrow(clusters_list[[1]][[nprime-1]]), ncol = ncol(data_wv))
    
    
    for(i in 1:nrow(quasi_U_statistics)){
      
      for(j in 1:ncol(quasi_U_statistics)){
        
        quasi_U_statistics[i,j] <- B_n(data_wv[clusters_list[[1]][[nprime-1]][i,],j],
                                       data_wv[clusters_list[[2]][[nprime-1]][i,],j],
                                       kernel)
        
      }
      
    }
    
    output[[nprime-1]] <- quasi_U_statistics
    
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
data_matrix_gen <- function(t,n,n1,model){
  
  data <- matrix(NA, nrow = n, ncol = t)
  
  gen_fun <- get(paste0("Model_", model, "_sim"))
  
  for(i in 1:n){
    
    if(i <= n1){
      
      data[i,] <- gen_fun(t)
      
    } else{
      
      data[i,] <- 1.25*gen_fun(t)
      
    }
    
  }
  
  return(data)
  
}

#' Find maximum value and index per column
#' 
#' @param X Data matrix
#' @return Matrix with [max_row, col_index, max_value]
max_per_col <- function(X){
  
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
max_per_col_list <- function(Bn_matrices,n){
  
  max_Bn_scales <- vector(mode = "list", length = floor(n/2)-1)
  
  for(nprime in 2:floor(n/2)){
    
    max_Bn_scales[[nprime-1]] <- max_per_col(Bn_matrices[[nprime-1]])
    
  }
  
  return(max_Bn_scales)
  
}

#' Bootstrap Bn statistic across scales
#' 
#' @param data_wv_list List containing wavelet variances for two groups
#' @return A vector of pooled bootstrap Bn estimates
Bootstrap_Bn_multiscale <- function(data_wv_list){
  
  data1 <- data_wv_list[[1]]
  data2 <- data_wv_list[[2]]
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n <- n1 + n2
  
  data <- rbind(data1, data2)
  
  indexes <- sample(n, replace = TRUE)
  
  data <- data[indexes,]
  
  data1_temp <- data[1:n1,]
  data2_temp <- data[(n1+1):n,]
  
  #return(list(data1_temp, data2_temp))
  
  output <- NULL
  
  for(j in 1:ncol(data1)){
    
    output <- c(output,B_n(data1_temp[,j],
                           data2_temp[,j],
                           kernel)) 
    
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
Bootstrap_Bn_multiscale_list <- function(data_wv, clusters_list, n, B){
  
  Bootstrap_Bn_versions <- vector(mode = "list", length = floor(n/2)-1)
  
  for(b in 1:B){
    
    for(nprime in 2:floor(n/2)){
      
      data_wv_list <- list(data_wv[clusters_list[[1]][[nprime-1]][1,],],
                           data_wv[clusters_list[[2]][[nprime-1]][1,],])
      
      Bootstrap_Bn_versions[[nprime-1]] <- rbind(Bootstrap_Bn_versions[[nprime-1]],
                                                 Bootstrap_Bn_multiscale(data_wv_list))
      
    }
    
  }
  
  return(Bootstrap_Bn_versions)
  
}

#' Estimate standard errors from bootstrap samples
#' 
#' @param Bootstrap_Bn_versions List of bootstrap Bn matrices
#' @param n Total elements
#' @return A list of standard error vectors
std.errors_calc <- function(Bootstrap_Bn_versions, n){
  
  std.errors <- vector(mode = "list", length = floor(n/2)-1)
  
  for(nprime in 2:floor(n/2)){
    
    std.errors[[nprime-1]] <- apply(Bootstrap_Bn_versions[[nprime-1]], 2, function(x){sd(x)})
    
  }
  
  return(std.errors)
  
}

#' Calculate p-values for observed Bn maxima
#' 
#' @param max_Bn_scales List of observed maxima
#' @param std.errors List of estimated standard errors
#' @param n Total elements
#' @return List of p-values
p_values_calc <- function(max_Bn_scales, std.errors, n){
  
  p_values <- vector(mode = "list", length = floor(n/2)-1)
  
  for(nprime in 2:floor(n/2)){
    
    for(j in 1:nrow(max_Bn_scales[[1]])){
      
      p_values[[nprime-1]] <- c(p_values[[nprime-1]], 1-pnorm(max_Bn_scales[[nprime-1]][j,3],sd = std.errors[[nprime-1]][j]))
      
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
best_clusters_fun <- function(p_values, t, max_Bn_scales, clusters_list){
  
  min_value <- Inf
  min_index <- NULL
  
  max_scale <- floor(log2(1+(t-1)/(8-1)))

  
  best_cluster <- vector(mode = "list", length = max_scale)
  
  for(scal in 1:max_scale){
    
    for (i in 1:length(p_values)){
      current_array <- p_values[[i]][1:scal]
      
      current_min <- min(current_array)
      current_min_index <- which(current_array == current_min, arr.ind = TRUE)
      
      if (current_min < min_value) {
        min_value <- current_min
        min_index <- c(i, current_min_index)
      }
      
    }
    
    best_cluster[[scal]] <- list(clusters_list[[1]][[min_index[1]]][max_Bn_scales[[min_index[[1]]]][min_index[[2]],1],],
                                 clusters_list[[2]][[min_index[1]]][max_Bn_scales[[min_index[[1]]]][min_index[[2]],1],])
    
    
  }
  
  
  return(best_cluster)
  
}

#' Convert cluster list to membership notation
#' 
#' @param cluster_list List of element indexes for each cluster
#' @param n Total elements
#' @return Vector of membership labels
cluster_notation_conv <- function(cluster_list, n){
  
  output <- rep(0,n)
  
  for(i in 1:length(cluster_list)){
    
    output[cluster_list[[i]]] <- i
    
  }
  
  return(output)
  
}


#' Run clustering analysis on real data
#' 
#' @param data Matrix of time series
#' @param B Bootstrap reps
#' @return A list of clustering results for ACF, COR, EUCL, PER and pooled Bn
cluster_fun <- function(data,B){
  
  n <- nrow(data)
  
  t <- ncol(data)
  
  data_wv <- t(apply(data, 1, wv_estimates))

  clusters_list <- possible_clusters(n)
  
  Bn_matrices <- Bn_matrices_list(data_wv, n, clusters_list)
  
  max_Bn_scales <- max_per_col_list(Bn_matrices,n)
  
  Bootstrap_Bn_versions <- Bootstrap_Bn_multiscale_list(data_wv, clusters_list, n, B)
  
  std.errors <- std.errors_calc(Bootstrap_Bn_versions, n)
  
  p_values <- p_values_calc(max_Bn_scales, std.errors, n)
  
  best_clusters <- best_clusters_fun(p_values, t, max_Bn_scales, clusters_list)
  
  best_clusters <- lapply(best_clusters, function(x){cluster_notation_conv(x, n)})
  
  #ACF
  
  dist_obj_ACF <- diss(data, "ACF")
  
  clust_ACF <- pam(dist_obj_ACF, k = 2, diss = TRUE)
  
  output_ACF <- clust_ACF$clustering
  
  #COR
  
  dist_obj_COR <- diss(data, "COR")
  
  clust_COR <- pam(dist_obj_COR, k = 2, diss = TRUE)
  
  output_COR <- clust_COR$clustering
  
  #EUCL
  
  dist_obj_EUCL <- diss(data, "EUCL")
  
  clust_EUCL <- pam(dist_obj_EUCL, k = 2, diss = TRUE)
  
  output_EUCL <- clust_EUCL$clustering
  
  #PER
  
  dist_obj_PER <- diss(data, "PER")
  
  clust_PER <- pam(dist_obj_PER, k = 2, diss = TRUE)
  
  output_PER <- clust_PER$clustering
  
  output <- c(list(output_ACF, output_COR, output_EUCL, output_PER), best_clusters)
  
  names(output) <- c("ACF", "COR", "EUCL", "PER", seq(1,length(output)-4)) 
  
  return(output)
  
}

################################################################################

#Arrowhead data

# Note: Data file is not included in the repository. See thesis for sources.
arrowhead <- unname(unlist(read.table(file.path(BASE_PATH, "Dados", "Arrowhead", "ArrowHead.txt"))))

plot(arrowhead, type = "l")

inflect <- function(x, threshold = 1){
  up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  a    <- cbind(x,up,down)
  list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
}

minima <- c(1,inflect(arrowhead, 30)$minima, 1506)
lenghts <- NULL

for(i in 1:24){
 
  lenghts <- c(lenghts, minima[i+1]-1-minima[i]+1)
   
}

min(lenghts)

arrowhead_list <- vector(mode = "list", length = 24)

for(i in 1:24){
  
  temp <- arrowhead[minima[i]:(minima[i+1]-1)]
  
  if(length(temp) >= 63){
    
    temp <- temp[1:62]
    
  }
  
  arrowhead_list[[i]] <- temp
  
}

arrowhead_list_diff <- lapply(arrowhead_list, diff)

arrowhead_diff <- matrix(unlist(arrowhead_list_diff), nrow = 24, byrow = TRUE)[c(1:10,13:22),]

plot(c(arrowhead_diff[1,],arrowhead_diff[2,],arrowhead_diff[3,],arrowhead_diff[4,],arrowhead_diff[5,],
       arrowhead_diff[6,],arrowhead_diff[7,],arrowhead_diff[8,],arrowhead_diff[9,],arrowhead_diff[10,]), 
     type = "l")
abline(v = seq(0.5,620.5, by = 62), lty = "dashed")

plot(c(arrowhead_diff[11,],arrowhead_diff[12,],arrowhead_diff[13,],arrowhead_diff[14,],arrowhead_diff[15,],
       arrowhead_diff[16,],arrowhead_diff[17,],arrowhead_diff[18,],arrowhead_diff[19,],arrowhead_diff[20,]), 
     type = "l")
abline(v = seq(0.5,620.5, by = 62), lty = "dashed")

set.seed(0)
results_arrowhead <- cluster_fun(arrowhead_diff, B)

results_arrowhead <- as.data.frame(matrix(unlist(results_arrowhead), ncol = 20, byrow = TRUE), row.names = names(results_arrowhead))

results_arrowhead[c(2,3,4),] <- 3 - results_arrowhead[c(2,3,4),]

################################################################################

#INMET data

# Note: Data file is not included in the repository. See thesis for sources.
INMET <- read.csv(file.path(BASE_PATH, "Dados", "Selecao", "dados_INMET_processados.csv"))

INMET <- t(INMET[seq(1,nrow(INMET), by = 12),-c(1,2)])

INMET_diff <- unname(t(apply(INMET, 1, diff)))

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_aplicacoes_part_1.RData"))
set.seed(0)
results_INMET <- cluster_fun(INMET_diff, B)

results_INMET <- as.data.frame(matrix(unlist(results_INMET), ncol = 10, byrow = TRUE), row.names = names(results_INMET))

results_INMET[-c(1,4),] <- 3 - results_INMET[-c(1,4),]

results_INMET <- 3 - results_INMET

save.image(file.path(WORKSPACE_DIR, "12_Clustering_Bn_aplicacoes_final.RData"))
