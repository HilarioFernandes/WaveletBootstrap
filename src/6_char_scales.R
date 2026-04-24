# =============================================================================
# 6_char_scales.R
# =============================================================================
# Purpose  : Characteristic scales study (scales where wavelet variance is maximised).
# Chapter  : Thesis Chapter 3
# Inputs   : None
# Outputs  : Characteristic scale estimates and CIs.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

# Set and create output directory for plots
OUTPUT_PATH <- file.path(BASE_PATH, "Plots/Plots_6")
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)

################################################################################

#characteristic scales

#' Find level with maximum wavelet variance
#' 
#' @param wv_est Vector of wavelet variance estimates
#' @return Index of the maximum value
max_char_scale <- function(wv_est){
  
  return(unname(which.max(wv_est)))
  
}

if (FALSE) {
# N <- 2048
# Y <- Model_A_sim(N)
# wv_est <- wv_estimates(Y)
# plot(log2(wv_est))
#  
# max_char_scale(wv_est)
}


#' Point estimate of largest local characteristic scale (level)
#' 
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return Quadratic interpolation of the characteristic scale
Char_scale_est <- function(wv_est, j_0){

  log_wv_est <- log2(wv_est)
  
  if(j_0 == 1 || j_0 == length(wv_est)){
    return(-1)
  }
  
  #we assume j_0 > 1
  
  x <- (j_0-1):(j_0+1)
  
  y <- log_wv_est[(j_0-1):(j_0+1)]
  
  quadratic <- solve(cbind(1, x, x^2), y)
  
  return(-quadratic[2]/(2*quadratic[3]))
  
}

if (FALSE) {
# N <- 2048
# Y <- Model_A_sim(N)
# wv_est <- wv_estimates(Y)
# plot(log2(wv_est))
# 
# Char_scale_est(wv_est,max_char_scale(wv_est))
}



#' beta_hat calculation for quadratic interpolation
#' 
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return Vector of beta coefficients
beta_hat_est <- function(wv_est, j_0){
  
  H <- rbind(c(-1/2,0,1/2),c(1,-2,1)) 
  
  return(as.vector(H %*% log2(wv_est[(j_0 - 1):(j_0 + 1)])))
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_A_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# wv_est <- wv_estimates(Y)
# 
# beta_hat_est(wv_est, 7)
}

#' Function that calculates the terms \hat{s}
#'
#' @param W List with all wavelet coefficients
#' @param N Time series length
#' @param L Array with the filter lengths
#' @return A list with biased autocovariance estimates for each scale
biased_est_acvs <- function(W, N, L){
  
  n_levels <- length(L)
  
  output <- vector(mode = "list", length = n_levels)
  
  for(j in 1:n_levels){
    
    temp <- NULL
    
    for(tau in 0:(N-L[j])){
      
      sum <- 0
      
      for(t in (L[j]-1):(N-tau-1)){
        
        sum <- sum + W[[j]][t]*W[[j]][t+tau]
        
      }
      
      temp <- c(temp, sum/(N-L[j]+1))
      
    }
    
    output[[j]] <- temp
    
  }
  
  return(output)
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_A_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# coeffs <- modwt(Y, n.levels = n_levels)[1:n_levels]
# 
# biased_est_acvs(coeffs, N, L)
}

#' Function that calculates the matrix \Sigma_1 from eq (8)
#' 
#' @param s_est List of biased acvs estimates
#' @param wv_est Vector of wavelet variance estimates
#' @param N Time series length
#' @param L Array of filter lengths
#' @param j_0 Target level index
#' @return The Sigma1 matrix estimate
Sigma1_est <- function(s_est, wv_est, N, L, j_0){
  
  temp <- matrix(NA, nrow = 3, ncol = 3)
  
  for(k_1 in 1:3){
    
    for(k_2 in k_1:3){
      
      sum <- 0  

      for(tau in 1:(N-L[j_0+(k_2-2)])){
        
        #the lag window chosen is the same as in the paper
        
        sum <- sum + s_est[[j_0+(k_1-2)]][tau]*s_est[[j_0+(k_2-2)]][tau]
        
      }
      
      sum <- (2*sum + wv_est[j_0+(k_1-2)]*wv_est[j_0+(k_2-2)])/2
      
      temp[k_2,k_1] <- 2*sum/(N-L[k_1]+1)
      
    }
    
  }
  
  return(temp)
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_C_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# j_0 <- 3
# 
# coeffs <- modwt(Y, n.levels = n_levels)[1:n_levels]
# 
# s_est <- biased_est_acvs(coeffs, N, L)
# 
# wv_est <- wv_estimates(Y)
# 
# Sigma1_est(s_est, wv_est, N, L, j_0)
}




#' Function that calculates \Sigma_2 from eq (9)
#' 
#' @param sigma1_est The Sigma1 matrix estimate
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return The Sigma2 matrix estimate
Sigma2_est <- function(sigma1_est, wv_est, j_0){
  
  temp <- matrix(NA, nrow = 3, ncol = 3)
  
  for(k_1 in 1:3){
    
    for(k_2 in k_1:3){
      
      sum <- sigma1_est[k_2, k_1]/(wv_est[k_1 + j_0 - 2]*wv_est[k_2 + j_0 - 2]*log(2)^2)
      
      sum <- sum + 2*(sigma1_est[k_1, k_1]*sigma1_est[k_2, k_2] + sigma1_est[k_2, k_1]^2)/(wv_est[k_1 + j_0 - 2]^2*wv_est[k_2 + j_0 - 2]^2*log(2)^2)
      
      
      
      temp[k_1,k_2] <- sum
      
      if(k_1 != k_2){
        
        temp[k_2,k_1] <- sum
        
      }
      
    }
    
  }
  
  return(temp)
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_C_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# j_0 <- 3
# 
# coeffs <- modwt(Y, n.levels = n_levels)[1:n_levels]
# 
# s_est <- biased_est_acvs(coeffs, N, L)
# 
# wv_est <- wv_estimates(Y)
# 
# sigma1_est <- Sigma1_est(s_est, wv_est, N, L, j_0)
# 
# Sigma2_est(sigma1_est, wv_est, j_0)
}


#' Function that calculates the var-cov matrix associated to \hat{\beta}
#' 
#' @param sigma2_est The Sigma2 matrix estimate
#' @return The variance-covariance matrix for beta coefficients
Var_beta_est <- function(sigma2_est){
  
  H <- rbind(c(-1/2,0,1/2),c(1,-2,1)) 
  
  return(H %*% sigma2_est %*% t(H))
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_C_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# j_0 <- 3
# 
# coeffs <- modwt(Y, n.levels = n_levels)[1:n_levels]
# 
# s_est <- biased_est_acvs(coeffs, N, L)
# 
# wv_est <- wv_estimates(Y)
# 
# sigma1_est <- Sigma1_est(s_est, wv_est, N, L, j_0)
# 
# sigma2_est <- Sigma2_est(sigma1_est, wv_est, j_0)
# 
# Var_beta_est(sigma2_est)
}


#' Function that calculates the variance of \hat{\kappa} from eq (10)
#' 
#' @param var_beta_est Variance-covariance matrix of beta
#' @param beta_hat Vector of beta coefficient estimates
#' @return Variance of the characteristic scale estimate
Var_kappa <- function(var_beta_est, beta_hat){
  
  sum <- var_beta_est[1,1]/beta_hat[2]^2 + (beta_hat[1]^2*var_beta_est[2,2])/beta_hat[2]^4
  
  sum <- sum + (var_beta_est[1,1]*var_beta_est[2,2] + 2*var_beta_est[1,2])/beta_hat[2]^4
  
  sum <- sum + (3*beta_hat[1]^2*var_beta_est[2,2]^2)/beta_hat[2]^6
  
  sum <- sum - (2*beta_hat[1]*var_beta_est[1,2])/beta_hat[2]^3
  
  return(sum)
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_A_sim(N)
# 
# n_levels <- floor(log2(1+(N-1)/(8-1)))
# 
# L <- (2^(1:n_levels)-1)*(8-1)+1
# 
# coeffs <- modwt(Y, n.levels = n_levels)[1:n_levels]
# 
# s_est <- biased_est_acvs(coeffs, N, L)
# 
# wv_est <- wv_estimates(Y)
# 
# j_0 <- max_char_scale(wv_est)
# 
# plot(log2(wv_est))
# 
# sigma1_est <- Sigma1_est(s_est, wv_est, N, L, j_0)
# 
# sigma2_est <- Sigma2_est(sigma1_est, wv_est, j_0)
# 
# var_beta_est <- Var_beta_est(sigma2_est)
# 
# beta_hat <- beta_hat_est(wv_est, max_char_scale(wv_est))
# 
# Var_kappa(var_beta_est, beta_hat)
}

#' Calculate confidence intervals for characteristic scales
#' 
#' @param data Time series data vector
#' @param alpha Significance level
#' @return A vector with [lower, upper] limits
CI_calc <- function(data, alpha){
  
  N <- length(data)
  
  n_levels <- floor(log2(1+(N-1)/(8-1)))
  
  L <- (2^(1:n_levels)-1)*(8-1)+1
  
  coeffs <- modwt(data, n.levels = n_levels)[1:n_levels]
  
  s_est <- biased_est_acvs(coeffs, N, L)
  
  wv_est <- wv_estimates(data)
  
  j_0 <- max_char_scale(wv_est)

  if(j_0 == 1 || j_0 == n_levels){
    
    return(c(-1,-1))
    
  }
  
  sigma1_est <- Sigma1_est(s_est, wv_est, N, L, j_0)
  
  sigma2_est <- Sigma2_est(sigma1_est, wv_est, j_0)
  
  var_beta_est <- Var_beta_est(sigma2_est)
  
  beta_hat <- beta_hat_est(wv_est, j_0)
  
  var_kappa <- Var_kappa(var_beta_est, beta_hat)
  
  char_scale_est <- Char_scale_est(wv_est, j_0)
  
  if(is.nan(sqrt(var_kappa))){
    
    return(c(-1,-1))
    
  }
  
  return(c(char_scale_est + qnorm(alpha/2)*sqrt(var_kappa),
           char_scale_est + qnorm(1 - alpha/2)*sqrt(var_kappa)))
  
}

if (FALSE) {
# N <- 2048
# 
# Y <- Model_A_sim(N)
# 
# CI_calc(Y, 0.05)
}

################################################################################

#study for CI

set.seed(45)

B <- 100
iterations <- 100

char_scale_est_A <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_B <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_C <- list("128" = NULL, "512" = NULL, "2048" = NULL)
char_scale_est_D <- list("128" = NULL, "512" = NULL, "2048" = NULL)

{
  
  A.wv_cs_CI <- list("128" = c(NULL,NULL),
                        "512" = c(NULL,NULL),
                        "2048" = c(NULL,NULL))
  B.wv_cs_CI <- list("128" = c(NULL,NULL),
                     "512" = c(NULL,NULL),
                     "2048" = c(NULL,NULL))
  C.wv_cs_CI <- list("128" = c(NULL,NULL),
                     "512" = c(NULL,NULL),
                     "2048" = c(NULL,NULL))
  D.wv_cs_CI <- list("128" = c(NULL,NULL),
                     "512" = c(NULL,NULL),
                     "2048" = c(NULL,NULL))
  
}

{
  
  A.wv_SB_2 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  A.wv_SB_4 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  A.wv_SB_8 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  B.wv_SB_2 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  B.wv_SB_4 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  B.wv_SB_8 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  C.wv_SB_2 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  C.wv_SB_4 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  C.wv_SB_8 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  D.wv_SB_2 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  D.wv_SB_4 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
  D.wv_SB_8 <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))
  
}




for(iter in 1:iterations){
  
  print(iter)
  
  for(i in 1:3){
    
    #print(c(iter,i))
    
    N <- 2^((2*i-1)+6)
    
    n_levels <- floor(log2(1+(N-1)/(8-1)))
    
    #simulating from each model
    YA <- Model_A_sim(N)
    YB <- Model_B_sim(N)
    YC <- Model_C_sim(N)
    YD <- Model_D_sim(N)
    
    #calculating the point estimates of the characteristic scales
    wv_est_A <- wv_estimates(YA)
    wv_est_B <- wv_estimates(YB)
    wv_est_C <- wv_estimates(YC)
    wv_est_D <- wv_estimates(YD)
    
    char_scale_est_A[[i]] <- c(char_scale_est_A[[i]], Char_scale_est(wv_est_A,max_char_scale(wv_est_A)))
    char_scale_est_B[[i]] <- c(char_scale_est_B[[i]], Char_scale_est(wv_est_B,max_char_scale(wv_est_B)))
    char_scale_est_C[[i]] <- c(char_scale_est_C[[i]], Char_scale_est(wv_est_C,max_char_scale(wv_est_C)))
    char_scale_est_D[[i]] <- c(char_scale_est_D[[i]], Char_scale_est(wv_est_D,max_char_scale(wv_est_D)))
    
    #calculating the confidence intervals via the normal approximation
    A.wv_cs_CI[[i]] <- rbind(A.wv_cs_CI[[i]], CI_calc(YA, 0.05)) 
    B.wv_cs_CI[[i]] <- rbind(B.wv_cs_CI[[i]], CI_calc(YB, 0.05))
    C.wv_cs_CI[[i]] <- rbind(C.wv_cs_CI[[i]], CI_calc(YC, 0.05))
    D.wv_cs_CI[[i]] <- rbind(D.wv_cs_CI[[i]], CI_calc(YD, 0.05))
    
    #bootstrap wavelet estimates
    
    A.wv_SB_2[[i]][[iter]] <- bootstrap_wavelet(YA, "bw", TRUE, "SB", function(N){1/(2*log2(N))}, B)
    
    A.wv_SB_4[[i]][[iter]] <- bootstrap_wavelet(YA, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    A.wv_SB_8[[i]][[iter]] <- bootstrap_wavelet(YA, "bw", TRUE, "SB", function(N){1/(8*log2(N))}, B)
    
    B.wv_SB_2[[i]][[iter]] <- bootstrap_wavelet(YB, "bw", TRUE, "SB", function(N){1/(2*log2(N))}, B)
    
    B.wv_SB_4[[i]][[iter]] <- bootstrap_wavelet(YB, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    B.wv_SB_8[[i]][[iter]] <- bootstrap_wavelet(YB, "bw", TRUE, "SB", function(N){1/(8*log2(N))}, B)
    
    C.wv_SB_2[[i]][[iter]] <- bootstrap_wavelet(YC, "bw", TRUE, "SB", function(N){1/(2*log2(N))}, B)
    
    C.wv_SB_4[[i]][[iter]] <- bootstrap_wavelet(YC, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    C.wv_SB_8[[i]][[iter]] <- bootstrap_wavelet(YC, "bw", TRUE, "SB", function(N){1/(8*log2(N))}, B)
    
    D.wv_SB_2[[i]][[iter]] <- bootstrap_wavelet(YD, "bw", TRUE, "SB", function(N){1/(2*log2(N))}, B)
    
    D.wv_SB_4[[i]][[iter]] <- bootstrap_wavelet(YD, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    D.wv_SB_8[[i]][[iter]] <- bootstrap_wavelet(YD, "bw", TRUE, "SB", function(N){1/(8*log2(N))}, B)
    
  }
  
}

################################################################################

true_char_scales <- function(list1){
  
  output <- NULL
  
  for(i in 1:3){
    
    if(length(which(list1[[i]] == -1)) > iterations/2){
      
      output <- c(output, -1)
      
    } else{
      
      temp <- list1[[i]]
      
      output <- c(output, median(temp[temp > 0]))
      
    }
    
  }
  
  return(output)
  
}

char_scale_A <- true_char_scales(char_scale_est_A)
char_scale_B <- true_char_scales(char_scale_est_B)
char_scale_C <- true_char_scales(char_scale_est_C)
char_scale_D <- true_char_scales(char_scale_est_D)



Char_scale_est_boot <- function(list1){
  
  temp <- list("128" = vector("list", length = iterations),
               "512" = vector("list", length = iterations),
               "2048" = vector("list", length = iterations))
  
  temp_fun <- function(x){Char_scale_est(x, max_char_scale(x))}
  
  for(i in 1:3){
    
    for(iter in 1:iterations){
      
      temp[[i]][[iter]] <- apply(list1[[i]][[iter]], 1, temp_fun)
      
    }
    
  }
  
  return(temp)

}

char_scale_est_A_SB_2_boot <- Char_scale_est_boot(A.wv_SB_2)
char_scale_est_A_SB_4_boot <- Char_scale_est_boot(A.wv_SB_4)
char_scale_est_A_SB_8_boot <- Char_scale_est_boot(A.wv_SB_8)

char_scale_est_B_SB_2_boot <- Char_scale_est_boot(B.wv_SB_2)
char_scale_est_B_SB_4_boot <- Char_scale_est_boot(B.wv_SB_4)
char_scale_est_B_SB_8_boot <- Char_scale_est_boot(B.wv_SB_8)

char_scale_est_C_SB_2_boot <- Char_scale_est_boot(C.wv_SB_2)
char_scale_est_C_SB_4_boot <- Char_scale_est_boot(C.wv_SB_4)
char_scale_est_C_SB_8_boot <- Char_scale_est_boot(C.wv_SB_8)

char_scale_est_D_SB_2_boot <- Char_scale_est_boot(D.wv_SB_2)
char_scale_est_D_SB_4_boot <- Char_scale_est_boot(D.wv_SB_4)
char_scale_est_D_SB_8_boot <- Char_scale_est_boot(D.wv_SB_8)

CI_calc_boot <- function(list1, alpha){
  
  output <- list(c(NULL,NULL),c(NULL,NULL),c(NULL,NULL))
  
  for(i in 1:3){
    
    for(iter in 1:iterations){
      
      if(length(which(list1[[i]][[iter]] == -1)) > B/2){
        
        output[[i]] <- rbind(output[[i]], c(-1,-1))
        
      } else{
        
        temp <- list1[[i]][[iter]]
        
        output[[i]] <- rbind(output[[i]], quantile(temp[temp > 0], c(alpha/2,1 - alpha/2)))
        
      }
     
    }
    
  }
  
  return(output)
  
}


A.SB_2_CI_boot <- CI_calc_boot(char_scale_est_A_SB_2_boot, 0.05)
A.SB_4_CI_boot <- CI_calc_boot(char_scale_est_A_SB_4_boot, 0.05)
A.SB_8_CI_boot <- CI_calc_boot(char_scale_est_A_SB_8_boot, 0.05)

B.SB_2_CI_boot <- CI_calc_boot(char_scale_est_B_SB_2_boot, 0.05)
B.SB_4_CI_boot <- CI_calc_boot(char_scale_est_B_SB_4_boot, 0.05)
B.SB_8_CI_boot <- CI_calc_boot(char_scale_est_B_SB_8_boot, 0.05)

C.SB_2_CI_boot <- CI_calc_boot(char_scale_est_C_SB_2_boot, 0.05)
C.SB_4_CI_boot <- CI_calc_boot(char_scale_est_C_SB_4_boot, 0.05)
C.SB_8_CI_boot <- CI_calc_boot(char_scale_est_C_SB_8_boot, 0.05)

D.SB_2_CI_boot <- CI_calc_boot(char_scale_est_D_SB_2_boot, 0.05)
D.SB_4_CI_boot <- CI_calc_boot(char_scale_est_D_SB_4_boot, 0.05)
D.SB_8_CI_boot <- CI_calc_boot(char_scale_est_D_SB_8_boot, 0.05)



CI_calc_boot_perc <- function(list1, array2){
  
  output <- NULL
  
  for(i in 1:3){
    
    temp_fun <- function(x){(x[1] <= array2[i]) & (x[2] >= array2[i])}
    
    output <- c(output, sum(apply(list1[[i]], 1, temp_fun))/iterations)
    
  }
  
  return(output)
  
}

A.SB_2_CI_boot_perc <- CI_calc_boot_perc(A.SB_2_CI_boot, char_scale_A)
A.SB_4_CI_boot_perc <- CI_calc_boot_perc(A.SB_4_CI_boot, char_scale_A)
A.SB_8_CI_boot_perc <- CI_calc_boot_perc(A.SB_8_CI_boot, char_scale_A)

B.SB_2_CI_boot_perc <- CI_calc_boot_perc(B.SB_2_CI_boot, char_scale_B)
B.SB_4_CI_boot_perc <- CI_calc_boot_perc(B.SB_4_CI_boot, char_scale_B)
B.SB_8_CI_boot_perc <- CI_calc_boot_perc(B.SB_8_CI_boot, char_scale_B)

C.SB_2_CI_boot_perc <- CI_calc_boot_perc(C.SB_2_CI_boot, char_scale_C)
C.SB_4_CI_boot_perc <- CI_calc_boot_perc(C.SB_4_CI_boot, char_scale_C)
C.SB_8_CI_boot_perc <- CI_calc_boot_perc(C.SB_8_CI_boot, char_scale_C)

D.SB_2_CI_boot_perc <- CI_calc_boot_perc(D.SB_2_CI_boot, char_scale_D)
D.SB_4_CI_boot_perc <- CI_calc_boot_perc(D.SB_4_CI_boot, char_scale_D)
D.SB_8_CI_boot_perc <- CI_calc_boot_perc(D.SB_8_CI_boot, char_scale_D)

A.wv_cs_CI_perc <- CI_calc_boot_perc(A.wv_cs_CI, char_scale_A)
B.wv_cs_CI_perc <- CI_calc_boot_perc(B.wv_cs_CI, char_scale_B)
C.wv_cs_CI_perc <- CI_calc_boot_perc(C.wv_cs_CI, char_scale_C)
D.wv_cs_CI_perc <- CI_calc_boot_perc(D.wv_cs_CI, char_scale_D)

results <- rbind(A.wv_cs_CI_perc,
              A.SB_2_CI_boot_perc,
              A.SB_4_CI_boot_perc,
              A.SB_8_CI_boot_perc,
              B.wv_cs_CI_perc,
              B.SB_2_CI_boot_perc,
              B.SB_4_CI_boot_perc,
              B.SB_8_CI_boot_perc,
              C.wv_cs_CI_perc,
              C.SB_2_CI_boot_perc,
              C.SB_4_CI_boot_perc,
              C.SB_8_CI_boot_perc,
              D.wv_cs_CI_perc,
              D.SB_2_CI_boot_perc,
              D.SB_4_CI_boot_perc,
              D.SB_8_CI_boot_perc)


for(i in 1:nrow(results)){
  
  str_temp <- ""
  
  if((i %% 4) == 1){
    
    str_temp <- c("(a)", "(b)", "(c)", "(d)")[1+(i %/% 4)]
    
  }
  
  cat(paste(c(str_temp, c("Apr. Normal", "$(2\\log_2(N))^{-1}$ ", "$(4\\log_2(N))^{-1}$ ",
                          "$(8\\log_2(N))^{-1}$ ")[((i-1) %% 4) + 1], results[i,]),
                          collapse = " & "))
  
  if(i < nrow(results)){
    
    cat(" \\\\")
    
  }
  
  if((i %% 4) == 0 && i != 16){
    
    cat(" \\hline")
    
  }
  
  cat("\n")
  
}

################################################################################

#Auxiliary plot

set.seed(1)

YA <- Model_A_sim(2048)
YB <- Model_B_sim(2048)
YC <- Model_C_sim(2048)
YD <- Model_D_sim(2048)

wv_est_A <- wv_estimates(YA)
wv_est_B <- wv_estimates(YB)
wv_est_C <- wv_estimates(YC)
wv_est_D <- wv_estimates(YD)


{
  png(file=file.path(OUTPUT_PATH, "Realizations.png"), width=1200, height=800, res = 150)

par(mfrow=c(2,2), mar = c(4.1, 4.1, 1.1, 2.1))

plot(log2(wv_est_A), xlab = "", ylab = "Var. de ondaletas (log)", pch = 19,
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
     lwd = 1.5, xaxt = 'n')
axis(1, labels = FALSE)

abline(v = char_scale_A[[3]], lty = 2,
       lwd = 1.5)

x <- 5:7
y <- log2(wv_est_A)[5:7]
quadratic <- solve(cbind(1, x, x^2), y)

quadratic_fun <- function(x){
  
  quadratic[[1]] + quadratic[[2]]*x + quadratic[[3]]*x^2
  
}

lines(seq(5,7,length.out = 50), quadratic_fun(seq(5,7,length.out = 50)),
      lwd = 1.5)

plot(log2(wv_est_B), xlab = "", ylab = "", pch = 19,
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
     lwd = 1.5, xaxt = 'n')
axis(1, labels = FALSE)

abline(v = char_scale_B[[3]], lty = 2,
       lwd = 1.5)

x <- 3:5
y <- log2(wv_est_B)[3:5]
quadratic <- solve(cbind(1, x, x^2), y)

lines(seq(3,5,length.out = 50), quadratic_fun(seq(3,5,length.out = 50)),
      lwd = 1.5)

plot(log2(wv_est_C), xlab = "Nível", ylab = "Var. de ondaletas (log)", pch = 19,
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
     lwd = 1.5)

abline(v = char_scale_C[[3]], lty = 2,
       lwd = 1.5)

x <- 4:6
y <- log2(wv_est_C[4:6])
quadratic <- solve(cbind(1, x, x^2), y)

lines(seq(4,6,length.out = 50), quadratic_fun(seq(4,6,length.out = 50)),
      lwd = 1.5)

plot(log2(wv_est_D), xlab = "Nível", ylab = "", pch = 19,
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
     lwd = 1.5)

dev.off()
}

################################################################################
