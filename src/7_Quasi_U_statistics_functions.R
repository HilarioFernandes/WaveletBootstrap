# =============================================================================
# 7_Quasi_U_statistics_functions.R
# =============================================================================
# Purpose  : Bn quasi U-statistic machinery for two-sample comparison.
# Chapter  : Chapter 4
# Inputs   : User-supplied symmetric kernel.
# Outputs  : Quasi U-statistic calculations and bootstrap/CI logic.
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================
#
# IMPORTANT: The calling environment must define 'kernel', e.g.:
# kernel <- function(x, y) (x - y)^2
# The kernel above was the one used in the simulation study.
#
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

#' Calculate U-statistic for a single group
#' 
#' @param data Group data vector
#' @param kernel Symmetric kernel function
#' @return The U-statistic estimate
U_ng <- function(data, kernel){
  
  n <- length(data)
  
  temp <- NULL
  
  for(i in 1:(n-1)){
    
    for(j in (i+1):n){
      
      temp <- c(temp, kernel(data[i], data[j]))
      
    }
    
  }
  
  return(mean(temp))
  
}

if (FALSE) {
# data <- rnorm(20)
# kernel <- function(x,y){
#   
#   return((x-y)^2)
#   
# }
# 
# U_ng(data, kernel)
}




#' Calculate cross-group U-statistic
#' 
#' @param data1 Group 1 data vector
#' @param data2 Group 2 data vector
#' @param kernel Symmetric kernel function
#' @return The cross-group U-statistic estimate
U_ngg <- function(data1,data2,kernel){
  
  n1 <- length(data1)
  n2 <- length(data2)
  
  temp <- NULL
  
  for(i in 1:n1){
    
    for(j in 1:n2){
      
      temp <- c(temp, kernel(data1[i], data2[j]))
      
    }
    
  }
  
  return(mean(temp))
  
  
}

if (FALSE) {
# data1 <- rnorm(20)
# data2 <- rnorm(20,1,1)
# 
# kernel <- function(x,y){
#   
#   return((x-y)^2)
#   
# }
# 
# U_ngg(data1, data2, kernel)
}









#' Calculate the Bn quasi U-statistic for two-sample comparison
#' 
#' @param data1 Group 1 data vector
#' @param data2 Group 2 data vector
#' @param kernel Symmetric kernel function
#' @return The Bn statistic estimate (n1*n2*(2*U_ngg - U_ng1 - U_ng2)/(n*(n-1)))
B_n <- function(data1, data2, kernel){
  
  n1 <- length(data1)
  n2 <- length(data2)
  n <- n1 + n2
  
  return( n1*n2*(2*U_ngg(data1, data2, kernel) - U_ng(data1, kernel) - U_ng(data2, kernel)  )/(n*(n-1)) )
  
}

if (FALSE) {
# data1 <- rnorm(20)
# data2 <- rnorm(25,1,1)
# 
# kernel <- function(x,y){
#   
#   return((x-y)^2)
#   
# }
# 
# B_n(data1, data2, kernel)
}



#' Bootstrap groups for Bn statistic
#' 
#' @param data_list List containing data1 and data2
#' @param discrimination Logical. If TRUE, resamples groups separately. 
#'        If FALSE, resamples from the pooled data (under H0).
#' @return A list with resampled data1_temp and data2_temp
bootstrap_Bn <- function(data_list, discrimination){
  
  data1 <- data_list[[1]]
  data2 <- data_list[[2]]

  n1 <- length(data1)
  n2 <- length(data2)
  n <- n1 + n2
    
  if(discrimination == TRUE){
    
    indexes1 <- sample(n1, replace = TRUE)
    indexes2 <- sample(n2, replace = TRUE)
    
    data1_temp <- data1[indexes1]
    data2_temp <- data2[indexes2]
    
  } else{
    
    data <- c(data1, data2)
    
    indexes <- sample(n, replace = TRUE)
    
    data <- data[indexes]
    
    data1_temp <- data[1:n1]
    data2_temp <- data[(n1+1):n]
    
  }
  
  return(list(data1_temp, data2_temp))
  
}

if (FALSE) {
# data1 <- runif(20, 0,1)
# data2 <- rnorm(25,-1,1)
# 
# bootstrap_Bn(list(data1, data2), TRUE)
# bootstrap_Bn(list(data1, data2), FALSE)
}


################################################################################  

#' Calculate quantile-based confidence interval for Bn
#' 
#' @param statistics_resampling Matrix of bootstrap statistics
#' @param alpha Significance level
#' @return Matrix of quantile CI limits
quantile_CI <- function(statistics_resampling, alpha){
  
  temp <- apply(statistics_resampling, 1, function(x){quantile(x, 1-alpha)})
  
  return(t(temp))
  
}

#' Calculate normal approximation confidence interval for Bn
#' 
#' @param statistics_resampling Matrix of resampling statistics
#' @param type Either "boot" or "jackknife"
#' @param alpha Significance level
#' @return Matrix of CI limits
#' @note For jackknife, variable 'n' must be in calling environment.
approx_CI <- function(statistics_resampling, type, alpha){
  
  if(type == "boot"){

    std.errors <- apply(statistics_resampling, 1, function(x){sd(x)})
    
  }else{
    
    std.errors <- apply(statistics_resampling, 1, function(x){sd(x)*(2*n-1)/(sqrt(2*n))})
    
  }
  
  return(t(sapply(std.errors, function(x){qnorm(1-alpha, mean = 0, sd = x)})))
  
}

#' Calculate coverage (rejection rate) of CI
#' 
#' @param statistics Vector of observed statistics
#' @param CI Matrix of CI limits
#' @return Rejection rate (1 means rejected, 0 not rejected)
coverage_CI <- function(statistics, CI){
  
  interval_check <- function(x){
    
    if(x[1] <= x[2]){
      return(0)
    } else{
      return(1)
    }
    
  }
  
  return(sum(apply(cbind(statistics, CI), 1, interval_check))/length(statistics))
  
}
