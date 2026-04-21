# =============================================================================
# 10_Comparison_Changepoint_functions.R
# =============================================================================
# Purpose  : Changepoint detection machinery (sequential and non-sequential).
# Chapter  : Thesis Chapter 3
# Inputs   : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Outputs  : Block index pairs indicating detected variance changes.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

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

#' Compare two blocks of data for variance change
#' 
#' @param X Data vector for block 1
#' @param Y Data vector for block 2
#' @param approximation Either "F" (F-test) or "boot_quant" (Bootstrap quantiles)
#' @param alpha_0 Significance level
#' @param jmax Maximum wavelet level to consider
#' @param B Number of bootstrap reps (if approximation == "boot_quant")
#' @return 1 if change detected, 0 otherwise
blocks_comparison <- function(X, Y, approximation, alpha_0, jmax, B){
  
  X.wv <- wv_estimates(X)
  Y.wv <- wv_estimates(Y)
  
  ratios <- Y.wv/X.wv
  
  X_Y <- c(X,Y)
  
  N1 <- length(X)
  N2 <- length(Y)
  N <- N1+N2
  
  alpha <- alpha_0/jmax
  
  if(approximation == "boot_quant"){
    
    ratios_boot <- NULL
    
    for(b in 1:B){
      
      boot_indexes <- block_boot(N, "SB", 1/(4*log2(N)))
      
      X_Y_boot <- X_Y[boot_indexes]
      X_boot.wv <- wv_estimates(X_Y_boot[1:N1])
      Y_boot.wv <- wv_estimates(X_Y_boot[(N1+1):(N)])
      
      ratios_boot <- rbind(ratios_boot, array(Y_boot.wv/X_boot.wv))
      
    } 
    
    quants <- t(apply(ratios_boot, 2, function(x){quantile(x,c(alpha/2,1-alpha/2))})[,1:jmax])
    
    rejection <- max(apply(cbind(ratios[1:jmax], quants), 1, interval_check))
    
  } else if(approximation == "F"){
    
    js <- 1:jmax
    
    L_js <- (2^(js) -1)*(8 - 1)+1
    
    N_js <- N - L_js + 1
    
    eta_js <- sapply(N_js/2^js, function(x){max(x,1)})
    
    quants <- t(sapply(eta_js, function(x){qf(c(alpha/2,1-alpha/2), x, x)}))
    
    rejection <- max(apply(cbind(ratios[1:jmax], quants), 1, interval_check))
    
  }
  
  return(rejection)
  
}


if (FALSE) {
# mult_factor <- 1
# 
# N <- 128
# X <- Model_B_sim(2*N)
# Y <- mult_factor*X[(N+1):(2*N)]
# X <- X[1:N]
# 
# jmax <- 3
# 
# alpha_0 <- 0.05
# 
# plot(c(X,Y), type = "l")
# 
# blocks_comparison(X, Y, "F", 0.05, 3, NA)
# blocks_comparison(X, Y, "boot_quant", 0.05, 3, 100)
}

################################################################################

#' Detect changepoints in a time series using block comparisons
#' 
#' @param data Time series data vector
#' @param block_size Length of each block
#' @param sequential Logical. If TRUE, compares current block to all previous blocks.
#'        If FALSE, only compares to the immediate previous block.
#' @param approximation Either "F" or "boot_quant"
#' @param alpha_0 Significance level
#' @param jmax Maximum wavelet level to consider
#' @param B Number of bootstrap reps
#' @return Matrix of block index pairs indicating where rejections occurred
changepointdetection <- function(data, block_size, sequential, approximation, alpha_0, jmax, B){
  
  amount_blocks <- as.integer(length(data)/block_size)
  
  changepoint_blocks <- NULL
  
  if(sequential == TRUE){
    
    for(index1 in 2:amount_blocks){
      
      for(index2 in (index1 - 1):1){
        
        X <- data[(block_size*(index1-1)+1):(index1*block_size)]
        Y <- data[(block_size*(index2-1)+1):(index2*block_size)]
        
        alpha <- alpha_0*(2^(-index2))/(1-2^(1-index1))
        
        rejection_temp <- blocks_comparison(Y,X,approximation,alpha, jmax, B) 
        
        if(rejection_temp == 1){
          
          changepoint_blocks <- rbind(changepoint_blocks, c(index2,index1))
          
          break
          
        }
        
      }
      
    }
    
  } else if (sequential == FALSE){
    
    for(index1 in 2:amount_blocks){
      
      index2 <- index1 - 1
      
      X <- data[(block_size*(index1-1)+1):(index1*block_size)]
      Y <- data[(block_size*(index2-1)+1):(index2*block_size)]
      
      rejection_temp <- blocks_comparison(Y, X, approximation, alpha_0, jmax, B) 
      
      if(rejection_temp == 1){
        
        changepoint_blocks <- rbind(changepoint_blocks, c(index2,index1))
        
      }
      
    }
    
  }
  
  return(changepoint_blocks)
  
}

if (FALSE) {
# #Example with two blocks of size 128
# 
# {
#   mult_factor <- 1
#   
#   N <- 128
#   X <- Model_B_sim(2*N)
#   Y <- mult_factor*X[(N+1):(2*N)]
#   X <- X[1:N]
#   
#   jmax <- 3
#   
#   alpha_0 <- 0.05
#   
#   plot(c(X,Y), type = "l")
# }
# 
# blocks_comparison(X, Y, "F", 0.05, 3, NA)
# changepointdetection(c(X,Y), 128, FALSE, "F", 0.05, 3, NA)
# changepointdetection(c(X,Y), 128, TRUE, "F", 0.05, 3, NA)
# 
# blocks_comparison(X, Y, "boot_quant", 0.05, 3, 100)
# changepointdetection(c(X,Y), 128, FALSE, "boot_quant", 0.05, 3, 100)
# changepointdetection(c(X,Y), 128, TRUE, "boot_quant", 0.05, 3, 100)
# 
# #Example with three blocks of size 128
# 
# {
#   N <- 128
#   X <- Model_B_sim(3*N)
#   Z <- X[(2*N+1):(3*N)]
#   Y <- X[(N+1):(2*N)]
#   X <- X[1:N]
#   
#   jmax <- 3
#   
#   alpha_0 <- 0.05
#   
#   plot(c(X,Y), type = "l")
# }
# 
# changepointdetection(c(X,Y,Z), 128, FALSE, "F", 0.05, 3, NA)
# changepointdetection(c(X,Y,Z), 128, TRUE, "F", 0.05, 3, NA)
# 
# changepointdetection(c(X,Y,Z), 128, FALSE, "boot_quant", 0.05, 3, 100)
# changepointdetection(c(X,Y,Z), 128, TRUE, "boot_quant", 0.05, 3, 100)
}


