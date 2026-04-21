# =============================================================================
# 9_Comparison_sim_aplicacoes.R
# =============================================================================
# Purpose  : Real-data comparison of feature sets (Arrowhead, INMET).
# Chapter  : Thesis Chapter 3
# Inputs   : ArrowHead.txt, dados_INMET_processados.csv (Data files not included).
# Outputs  : LaTeX tables for feature comparisons.
# Depends  : 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

B <- 100

alpha_0 <- 0.05

#' Check if ratio falls outside confidence interval
#' 
#' @param x Vector containing [ratio, lower_bound, upper_bound]
#' @return 1 if rejected (outside CI), 0 otherwise
interval_check <- function(x){
  
  if((x[1] <= x[2]) | (x[1] >= x[3])){
    
    return(1)
    
  } else{
    
    return(0)
    
  }
  
}

#' Convert binary array to range string
#' 
#' @param arr Binary vector
#' @return String representation of active indexes (e.g., "1-2, 4")
index_conversion <- function(arr) {
  # Find the indexes of 1s
  indexes <- which(arr == 1)
  
  # If there are no 1s, return NA
  if (length(indexes) == 0) {
    return(NA)
  }
  
  groups <- split(indexes, cumsum(c(1, diff(indexes) != 1)))
  
  result <- sapply(groups, function(g) {
    if (length(g) == 1) {
      return(as.character(g))
    } else {
      return(paste0(g[1], "-", g[length(g)]))
    }
  })
  
  result_string <- paste(result, collapse = ",")
  
  return(result_string)
}

#print(index_conversion(c(0,0,0,1))) # "4"
#print(index_conversion(c(1,1,0,0))) # "1-2"
#print(index_conversion(c(1,1,0,1))) # "1-2,4"
#print(index_conversion(c(1,1,1,1))) # "1-4"
#print(index_conversion(c(0,0,0,0))) # NA

#' Compare two time series features via wavelet variance ratios
#' 
#' @param X Data vector for group 1
#' @param Y Data vector for group 2
#' @param alpha Significance level
#' @param B Number of bootstrap reps
#' @return Matrix of rejections for F-test and Bootstrap (dep/indep)
comparison_fun <- function(X,Y, alpha, B){
  
  #we assume both have the same length
  N <- length(X)
  
  n_levels <- floor(log2(1+(N-1)/(8-1)))
  
  X.wv <- wv_estimates(X)
  Y.wv <- wv_estimates(Y)
  
  ratios <- Y.wv/X.wv
  
  X_Y <- c(X,Y)
  
  ratios_F <- NULL
  
  rejections_F <- NULL
  
  for(jmax in 1:n_levels){
    
    alpha <- alpha_0/jmax
    
    js <- 1:jmax
    
    L_js <- (2^(js) -1)*(8 - 1)+1
    
    N_js <- N - L_js + 1
    
    eta_js <- sapply(N_js/2^js, function(x){max(x,1)})
    
    quants_F <- t(sapply(eta_js, function(x){qf(c(alpha/2,1-alpha/2), x, x)}))
    rejections_F <- c(rejections_F,max(apply(cbind(ratios[1:jmax], quants_F), 1, interval_check)))
    
  }
  
  #bootstrap-based tests
  
  ratios_boot_indep <- NULL
  ratios_boot_dep <- NULL
  
  for(b in 1:B){
    
    #dep method
    
    boot_indexes <- block_boot(2*N, "SB", 1/(4*log2(2*N)))
    
    X_Y_boot <- X_Y[boot_indexes]
    X_boot.wv_dep <- wv_estimates(X_Y_boot[1:N])
    Y_boot.wv_dep <- wv_estimates(X_Y_boot[(N+1):(2*N)])
    
    #indep method
    
    boot_indexes_A <- block_boot(N, "SB", 1/(4*log2(N)))
    boot_indexes_B <- block_boot(N, "SB", 1/(4*log2(N)))
    
    X_boot <- X[boot_indexes_A]
    Y_boot <- Y[boot_indexes_B]
    
    X_boot.wv_indep <- wv_estimates(X_boot)
    Y_boot.wv_indep <- wv_estimates(Y_boot)
    
    #calculating
    
    ratios_boot_dep <- rbind(ratios_boot_dep,array(Y_boot.wv_dep/X_boot.wv_dep))
    ratios_boot_indep <- rbind(ratios_boot_indep,array(Y_boot.wv_indep/X_boot.wv_indep))
    
  }
  
  
  rejections_boot_dep <- NULL
  rejections_boot_indep <- NULL
  
  for(jmax in 1:n_levels){
    
    alpha <- alpha_0/jmax
    
    quants_boot_dep <- t(apply(ratios_boot_dep, 2, function(x){quantile(x,c(alpha/2,1-alpha/2))})[,1:jmax])
    rejections_boot_dep <- c(rejections_boot_dep, max(apply(cbind(ratios[1:jmax], quants_boot_dep), 1, interval_check)))
    
    quants_boot_indep <- t(apply(ratios_boot_indep, 2, function(x){quantile(x,c(alpha/2,1-alpha/2))})[,1:jmax])
    rejections_boot_indep <- c(rejections_boot_indep, max(apply(cbind(ratios[1:jmax], quants_boot_indep), 1, interval_check)))
    
  }
  
  return(cbind(rejections_F, rejections_boot_dep, rejections_boot_indep))
  
}

#' Perform multiple pairwise feature comparisons
#' 
#' @param data Matrix of time series datasets
#' @param alpha Significance level
#' @param B Number of bootstrap reps
#' @return A list of matrices (F, bootstrap_sem_sep, bootstrap_com_sep)
multiple_comparison_fun <- function(data, alpha, B){
  
  output1 <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  output2 <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  output3 <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  
  colnames(output1) <- colnames(data)
  rownames(output1) <- colnames(data)
  
  colnames(output2) <- colnames(data)
  rownames(output2) <- colnames(data)
  
  colnames(output3) <- colnames(data)
  rownames(output3) <- colnames(data)
  
  k <- 0
  
  for(i in 1:(ncol(data)-1)){
    
    for(j in (i+1):ncol(data)){
      
      message(round(2*k/(ncol(data)*(ncol(data)-1)),3),"\r",appendLF=FALSE)
      flush.console()
      Sys.sleep(0.1)
      
      k <- k+1
      
      X <- data[,i]
      Y <- data[,j]
      
      comparison_temp <- comparison_fun(X,Y,alpha,B)
      
      output1[i,j] <- index_conversion(comparison_temp[,1])
      output2[i,j] <- index_conversion(comparison_temp[,2])
      output3[i,j] <- index_conversion(comparison_temp[,3])
      
    }
    
  }
  
  return(list("F" = output1, "bootstrap_sem_sep" = output2, "bootstrap_com_sep" = output3))
  
}


################################################################################

#Interest rates data

# interest_rates <- read.csv("C:/Users/Hilar/Downloads/Interest Rates/Interest_rates.csv", row.names = 1)
# 
# plot(interest_rates[1:239,1], type = "l")
# lines(interest_rates[1:239,2], col = "red")
# lines(interest_rates[1:239,3], col = "blue")
# lines(interest_rates[1:239,4], col = "green")
# 
# interest_rates_diff <- apply(interest_rates[1:239,], 2, diff)
# 
# plot(interest_rates_diff[,1], type = "l")
# lines(interest_rates_diff[,2], col = "red")
# lines(interest_rates_diff[,3], col = "blue")
# lines(interest_rates_diff[,4], col = "green")
# 
# interest_rates_diff2 <- apply(interest_rates[1:239,], 2, function(x){diff(x, differences = 2)})
# 
# plot(interest_rates_diff2[,1], type = "l")
# lines(interest_rates_diff2[,2], col = "red")
# lines(interest_rates_diff2[,3], col = "blue")
# lines(interest_rates_diff2[,4], col = "green")
# 
# interest_rates_log <- apply(interest_rates[1:239,], 2, log)
# 
# plot(interest_rates_log[1:239,1], type = "l")
# lines(interest_rates_log[1:239,2], col = "red")
# lines(interest_rates_log[1:239,3], col = "blue")
# lines(interest_rates_log[1:239,4], col = "green")
# 
# 
# X <- interest_rates_diff[,1]
# Y <- interest_rates_diff[,3]
# 
# comparison_fun(X,Y,0.05,100)
# 
# set.seed(0)
# results_interest_rates <- multiple_comparison_fun(interest_rates_diff, 0.05, 100)
# results_interest_rates_nodiff <- multiple_comparison_fun(interest_rates[1:239,], 0.05, 100)
# 
# results_interest_rates_log <- multiple_comparison_fun(interest_rates_log[1:239,], 0.05, 100)

################################################################################

#Arrowhead

# Note: Data file is not included in the repository. See thesis for sources.
arrowhead <- unname(unlist(read.table(file.path(BASE_PATH, "Arrowhead", "ArrowHead.txt"))))

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


sample_arrowhead <- cbind(arrowhead_list_diff[[1]], arrowhead_list_diff[[2]],
                          arrowhead_list_diff[[3]], arrowhead_list_diff[[4]],
                          arrowhead_list_diff[[5]],arrowhead_list_diff[[13]], arrowhead_list_diff[[14]],
                          arrowhead_list_diff[[15]], arrowhead_list_diff[[16]],
                          arrowhead_list_diff[[17]])

set.seed(0)
results_arrowhead <- multiple_comparison_fun(sample_arrowhead, 0.05, 100)

################################################################################

#INMET data

# Note: Data file is not included in the repository. See thesis for sources.
INMET <- read.csv(file.path(BASE_PATH, "Dados", "Selecao", "dados_INMET_processados.csv"))

INMET_selec <- INMET[seq(1,8760, by = 12),] 

plot(INMET_selec[,3], type = "l")
lines(INMET_selec[,8], col = "blue")

INMET_selec_diff <- apply(INMET_selec[-(1:2)], 2, diff)

plot(INMET_selec_diff[,1], type = "l")
lines(INMET_selec_diff[,5], col = "blue")

comparison_fun(INMET_selec_diff[,1],
               INMET_selec_diff[,6],
               0.05,100)

set.seed(0)
results_INMET <- multiple_comparison_fun(INMET_selec_diff, 0.05, 100)


