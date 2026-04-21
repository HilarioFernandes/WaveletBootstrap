# =============================================================================
# 5_CI_wv.R
# =============================================================================
# Purpose  : Simulation study of bootstrap CI coverage for wavelet variance.
# Chapter  : Thesis Chapter 3
# Inputs   : None
# Outputs  : LaTeX tables for CI coverage.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))################################################################################

#study for CI

set.seed(44)

B <- 100
iterations <- 100

A.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
B.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
C.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
D.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

{
  
A.wv_waveslim_gaussian <- list("128" = vector("list", length = iterations),
             "512" = vector("list", length = iterations),
             "2048" = vector("list", length = iterations))
B.wv_waveslim_gaussian <- list("128" = vector("list", length = iterations),
             "512" = vector("list", length = iterations),
             "2048" = vector("list", length = iterations))
C.wv_waveslim_gaussian <- list("128" = vector("list", length = iterations),
             "512" = vector("list", length = iterations),
             "2048" = vector("list", length = iterations))
D.wv_waveslim_gaussian <- list("128" = vector("list", length = iterations),
             "512" = vector("list", length = iterations),
             "2048" = vector("list", length = iterations))

}

{
  
  A.wv_waveslim_nongaussian <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  B.wv_waveslim_nongaussian <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  C.wv_waveslim_nongaussian <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  D.wv_waveslim_nongaussian <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  
}

{
  
  A.wv_waveslim_eta3 <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  B.wv_waveslim_eta3 <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  C.wv_waveslim_eta3 <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  D.wv_waveslim_eta3 <- list("128" = vector("list", length = iterations),
                        "512" = vector("list", length = iterations),
                        "2048" = vector("list", length = iterations))
  
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
  
  message(iter/iterations,"\r",appendLF=FALSE)
  flush.console()
  Sys.sleep(0.1)
  
  for(i in 1:3){
    
    #print(c(iter,i))
    
    N <- 2^((2*i-1)+6)
    
    n_levels <- floor(log2(1+(N-1)/(8-1)))
    
    #simulating from each model
    YA <- Model_A_sim(N)
    YB <- Model_B_sim(N)
    YC <- Model_C_sim(N)
    YD <- Model_D_sim(N)
    
    #calculating the point estimates of the wavelet variances
    A.wv[[i]] <- rbind(A.wv[[i]], wv_estimates(YA))
    B.wv[[i]] <- rbind(B.wv[[i]], wv_estimates(YB))
    C.wv[[i]] <- rbind(C.wv[[i]], wv_estimates(YC))
    D.wv[[i]] <- rbind(D.wv[[i]], wv_estimates(YD))
    
    #calculating the confidence intervals via waveslim
    A.wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(YA, n.levels = n_levels),
                                                p = 0.025, type = "gaussian")[1:n_levels,]
    B.wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(YB, n.levels = n_levels),
                                                p = 0.025, type = "gaussian")[1:n_levels,]
    C.wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(YC, n.levels = n_levels),
                                                p = 0.025, type = "gaussian")[1:n_levels,]
    D.wv_waveslim_gaussian[[i]][[iter]] <- wave.variance(modwt(YD, n.levels = n_levels),
                                                p = 0.025, type = "gaussian")[1:n_levels,]
    
    A.wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(YA, n.levels = n_levels),
                                                         p = 0.025, type = "nongaussian")[1:n_levels,]
    B.wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(YB, n.levels = n_levels),
                                                         p = 0.025, type = "nongaussian")[1:n_levels,]
    C.wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(YC, n.levels = n_levels),
                                                         p = 0.025, type = "nongaussian")[1:n_levels,]
    D.wv_waveslim_nongaussian[[i]][[iter]] <- wave.variance(modwt(YD, n.levels = n_levels),
                                                         p = 0.025, type = "nongaussian")[1:n_levels,]
    
    
    A.wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(YA, n.levels = n_levels),
                                                         p = 0.025, type = "eta3")[1:n_levels,]
    B.wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(YB, n.levels = n_levels),
                                                         p = 0.025, type = "eta3")[1:n_levels,]
    C.wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(YC, n.levels = n_levels),
                                                         p = 0.025, type = "eta3")[1:n_levels,]
    D.wv_waveslim_eta3[[i]][[iter]] <- wave.variance(modwt(YD, n.levels = n_levels),
                                                         p = 0.025, type = "eta3")[1:n_levels,]
    
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

#confidence intervals

# CI_calculation_boot <- function(list1, alpha){
#   
#   lower <- list("128" = NULL, "512" = NULL, "2048" = NULL)
#   upper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
#   
#   for(iter in 1:iterations){
# 
#     for(i in 1:3){
#       
#       N <- 2^((2*i-1)+6)
#       
#       n_levels <- floor(log2(1+(N-1)/(8-1)))
#       
#       lower[[i]] <- rbind(lower[[i]], apply(list1[[i]][[iter]], 2, function(x){quantile(x, probs = alpha/2)}))
#       
#       upper[[i]] <- rbind(upper[[i]], apply(list1[[i]][[iter]], 2, function(x){quantile(x, probs = 1 - (alpha/2))}))
#       
#     }
#     
#   }
#   
#   return(list(lower, upper))
#   
# }

#' Calculate bootstrap confidence intervals
#' 
#' @param list1 List of bootstrap variances estimates
#' @param list2 List of point estimates for wavelet variances
#' @param alpha Significance level
#' @return A list containing lower and upper limits matrices
CI_calculation_boot <- function(list1, list2, alpha){
  
  lower <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  upper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  
  for(iter in 1:iterations){
    
    for(i in 1:3){
      
      N <- 2^((2*i-1)+6)
      
      n_levels <- floor(log2(1+(N-1)/(8-1)))
      
      #bootstrap variances estimates and point estimates for wavelet variances
      
      vars_boot <- apply(list1[[i]][[iter]], 2, function(x){var(x)})
      
      mean_estimates <- list2[[i]][iter,]
      
      #assuming that the wavelet variance estimators are proportional to a chi-squared
      #distribution, we may estimate these two parameters via
      
      a <- vars_boot/(2*mean_estimates)
      
      eta <- (2*mean_estimates^2)/vars_boot
      
      lower[[i]] <- rbind(lower[[i]], a*qchisq(rep(alpha/2, n_levels), df= eta))
      
      upper[[i]] <- rbind(upper[[i]], a*qchisq(rep(1 - (alpha/2), n_levels), df= eta))
    }
    
  }
  
  return(list(lower, upper))
  
}

A.wv_SB_2_CI <- CI_calculation_boot(A.wv_SB_2, A.wv, 0.05)
A.wv_SB_4_CI <- CI_calculation_boot(A.wv_SB_4, A.wv, 0.05)
A.wv_SB_8_CI <- CI_calculation_boot(A.wv_SB_8, A.wv, 0.05)

B.wv_SB_2_CI <- CI_calculation_boot(B.wv_SB_2, B.wv, 0.05)
B.wv_SB_4_CI <- CI_calculation_boot(B.wv_SB_4, B.wv, 0.05)
B.wv_SB_8_CI <- CI_calculation_boot(B.wv_SB_8, B.wv, 0.05)

C.wv_SB_2_CI <- CI_calculation_boot(C.wv_SB_2, C.wv, 0.05)
C.wv_SB_4_CI <- CI_calculation_boot(C.wv_SB_4, C.wv, 0.05)
C.wv_SB_8_CI <- CI_calculation_boot(C.wv_SB_8, C.wv, 0.05)

D.wv_SB_2_CI <- CI_calculation_boot(D.wv_SB_2, D.wv, 0.05)
D.wv_SB_4_CI <- CI_calculation_boot(D.wv_SB_4, D.wv, 0.05)
D.wv_SB_8_CI <- CI_calculation_boot(D.wv_SB_8, D.wv, 0.05)




#' Calculate true wavelet variances
#' 
#' @param list1 Nested list of point estimates for wavelet variances
#' @return A vector of true wavelet variances
true_wv <- function(list1){
  
  temp <- apply(list1[[3]], 2, function(x){mean(x)})
  
  return(temp)

}

A.wv_mean <- true_wv(A.wv)
B.wv_mean <- true_wv(B.wv)
C.wv_mean <- true_wv(C.wv)
D.wv_mean <- true_wv(D.wv)

################################################################################

#coverages

#' Check CI coverage for bootstrap methods
#' 
#' @param list1 List of logical coverage matrices
#' @param array2 Array of true values
#' @return List of coverage true/false
coverages_boot <- function(list1, array2){
  
  output <- list(matrix(-1,nrow = iterations, ncol = 4),
                 matrix(-1,nrow = iterations, ncol = 6),
                 matrix(-1,nrow = iterations, ncol = 8))
  
  for(i in 1:3){
    
    N <- 2^((2*i-1)+6)
    
    n_levels <- floor(log2(1+(N-1)/(8-1)))
    
    for(iter in 1:iterations){
      
      for(j in 1:n_levels){
        
        if(array2[j] >= list1[[1]][[i]][iter,j] & array2[j] <= list1[[2]][[i]][iter,j]){
          
          output[[i]][iter,j] <- TRUE
          
        } else if(array2[j] < list1[[1]][[i]][iter,j] | array2[j] > list1[[2]][[i]][iter,j]){
          
          output[[i]][iter,j] <- FALSE
          
        }
        
      }
      
    }
    
  }
  
  return(output)
  
}

A.wv_SB_2_CI_coverage <- coverages_boot(A.wv_SB_2_CI, A.wv_mean)
A.wv_SB_4_CI_coverage <- coverages_boot(A.wv_SB_4_CI, A.wv_mean)
A.wv_SB_8_CI_coverage <- coverages_boot(A.wv_SB_8_CI, A.wv_mean)

B.wv_SB_2_CI_coverage <- coverages_boot(B.wv_SB_2_CI, B.wv_mean)
B.wv_SB_4_CI_coverage <- coverages_boot(B.wv_SB_4_CI, B.wv_mean)
B.wv_SB_8_CI_coverage <- coverages_boot(B.wv_SB_8_CI, B.wv_mean)

C.wv_SB_2_CI_coverage <- coverages_boot(C.wv_SB_2_CI, C.wv_mean)
C.wv_SB_4_CI_coverage <- coverages_boot(C.wv_SB_4_CI, C.wv_mean)
C.wv_SB_8_CI_coverage <- coverages_boot(C.wv_SB_8_CI, C.wv_mean)

D.wv_SB_2_CI_coverage <- coverages_boot(D.wv_SB_2_CI, D.wv_mean)
D.wv_SB_4_CI_coverage <- coverages_boot(D.wv_SB_4_CI, D.wv_mean)
D.wv_SB_8_CI_coverage <- coverages_boot(D.wv_SB_8_CI, D.wv_mean)

#' Check CI coverage using waveslim
#' 
#' @param list1 List of waveslim CI matrices
#' @param array2 Array of true values
#' @return List of coverage true/false
coverages_waveslim <- function(list1, array2){
  
  output <- list(matrix(-1,nrow = iterations, ncol = 4),
                 matrix(-1,nrow = iterations, ncol = 6),
                 matrix(-1,nrow = iterations, ncol = 8))
  
  for(i in 1:3){
    
    N <- 2^((2*i-1)+6)
    
    n_levels <- floor(log2(1+(N-1)/(8-1)))
    
    for(iter in 1:iterations){
      
      for(j in 1:n_levels){
        
        if(array2[j] >= list1[[i]][[iter]][j,2] & array2[j] <= list1[[i]][[iter]][j,3]){
          
          output[[i]][iter,j] <- TRUE
          
        } else if(array2[j] < list1[[i]][[iter]][j,2] | array2[j] > list1[[i]][[iter]][j,2]){
          
          output[[i]][iter,j] <- FALSE
          
        }
        
      }
      
    }
    
  }
  
  return(output)
  
}

A.wv_waveslim_gaussian_coverage <- coverages_waveslim(A.wv_waveslim_gaussian, A.wv_mean)
B.wv_waveslim_gaussian_coverage <- coverages_waveslim(B.wv_waveslim_gaussian, B.wv_mean)
C.wv_waveslim_gaussian_coverage <- coverages_waveslim(C.wv_waveslim_gaussian, C.wv_mean)
D.wv_waveslim_gaussian_coverage <- coverages_waveslim(D.wv_waveslim_gaussian, D.wv_mean)

A.wv_waveslim_nongaussian_coverage <- coverages_waveslim(A.wv_waveslim_nongaussian, A.wv_mean)
B.wv_waveslim_nongaussian_coverage <- coverages_waveslim(B.wv_waveslim_nongaussian, B.wv_mean)
C.wv_waveslim_nongaussian_coverage <- coverages_waveslim(C.wv_waveslim_nongaussian, C.wv_mean)
D.wv_waveslim_nongaussian_coverage <- coverages_waveslim(D.wv_waveslim_nongaussian, D.wv_mean)

A.wv_waveslim_eta3_coverage <- coverages_waveslim(A.wv_waveslim_eta3, A.wv_mean)
B.wv_waveslim_eta3_coverage <- coverages_waveslim(B.wv_waveslim_eta3, B.wv_mean)
C.wv_waveslim_eta3_coverage <- coverages_waveslim(C.wv_waveslim_eta3, C.wv_mean)
D.wv_waveslim_eta3_coverage <- coverages_waveslim(D.wv_waveslim_eta3, D.wv_mean)

################################################################################

#coverages (%)

#' Calculate coverage percentage
#' 
#' @param list1 List of boolean coverages matrices
#' @return The coverage percentage vector
coverages_percentual <- function(list1){
  
  temp <- list(NULL,NULL,NULL)
  
  for(i in 1:3){
    
    temp[[i]] <- apply(list1[[i]], 2, function(x){sum(x)/length(x)})
    
  }
  
  return(temp)
  
}

A.wv_SB_2_CI_coverage_percentual <- coverages_percentual(A.wv_SB_2_CI_coverage)
A.wv_SB_4_CI_coverage_percentual <- coverages_percentual(A.wv_SB_4_CI_coverage)
A.wv_SB_8_CI_coverage_percentual <- coverages_percentual(A.wv_SB_8_CI_coverage)

B.wv_SB_2_CI_coverage_percentual <- coverages_percentual(B.wv_SB_2_CI_coverage)
B.wv_SB_4_CI_coverage_percentual <- coverages_percentual(B.wv_SB_4_CI_coverage)
B.wv_SB_8_CI_coverage_percentual <- coverages_percentual(B.wv_SB_8_CI_coverage)

C.wv_SB_2_CI_coverage_percentual <- coverages_percentual(C.wv_SB_2_CI_coverage)
C.wv_SB_4_CI_coverage_percentual <- coverages_percentual(C.wv_SB_4_CI_coverage)
C.wv_SB_8_CI_coverage_percentual <- coverages_percentual(C.wv_SB_8_CI_coverage)

D.wv_SB_2_CI_coverage_percentual <- coverages_percentual(D.wv_SB_2_CI_coverage)
D.wv_SB_4_CI_coverage_percentual <- coverages_percentual(D.wv_SB_4_CI_coverage)
D.wv_SB_8_CI_coverage_percentual <- coverages_percentual(D.wv_SB_8_CI_coverage)

A.wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(A.wv_waveslim_gaussian_coverage)
B.wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(B.wv_waveslim_gaussian_coverage)
C.wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(C.wv_waveslim_gaussian_coverage)
D.wv_waveslim_gaussian_coverage_percentual <- coverages_percentual(D.wv_waveslim_gaussian_coverage)

A.wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(A.wv_waveslim_nongaussian_coverage)
B.wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(B.wv_waveslim_nongaussian_coverage)
C.wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(C.wv_waveslim_nongaussian_coverage)
D.wv_waveslim_nongaussian_coverage_percentual <- coverages_percentual(D.wv_waveslim_nongaussian_coverage)

A.wv_waveslim_eta3_coverage_percentual <- coverages_percentual(A.wv_waveslim_eta3_coverage)
B.wv_waveslim_eta3_coverage_percentual <- coverages_percentual(B.wv_waveslim_eta3_coverage)
C.wv_waveslim_eta3_coverage_percentual <- coverages_percentual(C.wv_waveslim_eta3_coverage)
D.wv_waveslim_eta3_coverage_percentual <- coverages_percentual(D.wv_waveslim_eta3_coverage)

## ---- Tables and LateX Generation ----

#' Format simulation results to LaTeX rules
#' 
#' @param results Aggregated matrix of percentages
#' @return Latex table row printed to console
latex_fun <- function(results){
  
  for(i in 1:nrow(results)){
    
    str_temp <- ""
    
    if((i+1) %% 3 == 0){
      
      str_temp <- c("Gaussian ", "$\\hat{\\eta}_3$ ", "Multitaper ", "$(2\\log_2(N))^{-1}$ ",
                    "$(4\\log_2(N))^{-1}$ ", "$(8\\log_2(N))^{-1}$ ")[(i+1) %/% 3]
      
    }
    
    cat(paste(c(str_temp, c(128,512,2048)[((i-1) %% 3) + 1], results[i,]), collapse = " & "))
    
    if(i < nrow(results)){
      
      cat(" \\\\")
      
    }
    
    
    if((i %% 3) == 0 && i != 18){
      
      cat("\\hline")
      
    }
    
    cat("\n")
    
  }
  
  
}

resultsA <- rbind(c(A.wv_waveslim_gaussian_coverage_percentual[[1]], rep("",4)),
                  c(A.wv_waveslim_gaussian_coverage_percentual[[2]], rep("",2)),
                  A.wv_waveslim_gaussian_coverage_percentual[[3]],
                  c(A.wv_waveslim_eta3_coverage_percentual[[1]], rep("",4)),
                  c(A.wv_waveslim_eta3_coverage_percentual[[2]], rep("",2)),
                  A.wv_waveslim_eta3_coverage_percentual[[3]],
              c(A.wv_waveslim_nongaussian_coverage_percentual[[1]], rep("",4)),
              c(A.wv_waveslim_nongaussian_coverage_percentual[[2]], rep("",2)),
              A.wv_waveslim_nongaussian_coverage_percentual[[3]],
              c(A.wv_SB_2_CI_coverage_percentual[[1]], rep("",4)),
              c(A.wv_SB_2_CI_coverage_percentual[[2]], rep("",2)),
              A.wv_SB_2_CI_coverage_percentual[[3]],
              c(A.wv_SB_4_CI_coverage_percentual[[1]], rep("",4)),
              c(A.wv_SB_4_CI_coverage_percentual[[2]], rep("",2)),
              A.wv_SB_4_CI_coverage_percentual[[3]],
              c(A.wv_SB_8_CI_coverage_percentual[[1]], rep("",4)),
              c(A.wv_SB_8_CI_coverage_percentual[[2]], rep("",2)),
              A.wv_SB_8_CI_coverage_percentual[[3]])

resultsB <- rbind(c(B.wv_waveslim_gaussian_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_waveslim_gaussian_coverage_percentual[[2]], rep("",2)),
                  B.wv_waveslim_gaussian_coverage_percentual[[3]],
                  c(B.wv_waveslim_eta3_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_waveslim_eta3_coverage_percentual[[2]], rep("",2)),
                  B.wv_waveslim_eta3_coverage_percentual[[3]],
                  c(B.wv_waveslim_nongaussian_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_waveslim_nongaussian_coverage_percentual[[2]], rep("",2)),
                  B.wv_waveslim_nongaussian_coverage_percentual[[3]],
                  c(B.wv_SB_2_CI_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_SB_2_CI_coverage_percentual[[2]], rep("",2)),
                  B.wv_SB_2_CI_coverage_percentual[[3]],
                  c(B.wv_SB_4_CI_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_SB_4_CI_coverage_percentual[[2]], rep("",2)),
                  B.wv_SB_4_CI_coverage_percentual[[3]],
                  c(B.wv_SB_8_CI_coverage_percentual[[1]], rep("",4)),
                  c(B.wv_SB_8_CI_coverage_percentual[[2]], rep("",2)),
                  B.wv_SB_8_CI_coverage_percentual[[3]])

resultsC <- rbind(c(C.wv_waveslim_gaussian_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_waveslim_gaussian_coverage_percentual[[2]], rep("",2)),
                  C.wv_waveslim_gaussian_coverage_percentual[[3]],
                  c(C.wv_waveslim_eta3_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_waveslim_eta3_coverage_percentual[[2]], rep("",2)),
                  C.wv_waveslim_eta3_coverage_percentual[[3]],
                  c(C.wv_waveslim_nongaussian_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_waveslim_nongaussian_coverage_percentual[[2]], rep("",2)),
                  C.wv_waveslim_nongaussian_coverage_percentual[[3]],
                  c(C.wv_SB_2_CI_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_SB_2_CI_coverage_percentual[[2]], rep("",2)),
                  C.wv_SB_2_CI_coverage_percentual[[3]],
                  c(C.wv_SB_4_CI_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_SB_4_CI_coverage_percentual[[2]], rep("",2)),
                  C.wv_SB_4_CI_coverage_percentual[[3]],
                  c(C.wv_SB_8_CI_coverage_percentual[[1]], rep("",4)),
                  c(C.wv_SB_8_CI_coverage_percentual[[2]], rep("",2)),
                  C.wv_SB_8_CI_coverage_percentual[[3]])

resultsD <- rbind(c(D.wv_waveslim_gaussian_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_waveslim_gaussian_coverage_percentual[[2]], rep("",2)),
                  D.wv_waveslim_gaussian_coverage_percentual[[3]],
                  c(D.wv_waveslim_eta3_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_waveslim_eta3_coverage_percentual[[2]], rep("",2)),
                  D.wv_waveslim_eta3_coverage_percentual[[3]],
                  c(D.wv_waveslim_nongaussian_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_waveslim_nongaussian_coverage_percentual[[2]], rep("",2)),
                  D.wv_waveslim_nongaussian_coverage_percentual[[3]],
                  c(D.wv_SB_2_CI_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_SB_2_CI_coverage_percentual[[2]], rep("",2)),
                  D.wv_SB_2_CI_coverage_percentual[[3]],
                  c(D.wv_SB_4_CI_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_SB_4_CI_coverage_percentual[[2]], rep("",2)),
                  D.wv_SB_4_CI_coverage_percentual[[3]],
                  c(D.wv_SB_8_CI_coverage_percentual[[1]], rep("",4)),
                  c(D.wv_SB_8_CI_coverage_percentual[[2]], rep("",2)),
                  D.wv_SB_8_CI_coverage_percentual[[3]])

latex_fun(resultsA)
latex_fun(resultsB)
latex_fun(resultsC)
latex_fun(resultsD)

################################################################################


#' Detailed LaTeX formatting wrapper function
#' 
#' @param resultsA Results from Model A
#' @param resultsB Results from Model B
#' @param resultsC Results from Model C
#' @param resultsD Results from Model D
#' @return Latex table row median values printed to console
latex_fun2 <- function(resultsA, resultsB, resultsC, resultsD){
  
  for(i in 1:nrow(resultsA)){
  
    if(i %% 3 == 1){
      
      val1A <- median(as.numeric(resultsA[i,1:2]), na.rm = TRUE)
      val2A <- median(as.numeric(resultsA[i,3:4]), na.rm = TRUE)
      
      val1B <- median(as.numeric(resultsB[i,1:2]), na.rm = TRUE)
      val2B <- median(as.numeric(resultsB[i,3:4]), na.rm = TRUE)
      
      val1C <- median(as.numeric(resultsC[i,1:2]), na.rm = TRUE)
      val2C <- median(as.numeric(resultsC[i,3:4]), na.rm = TRUE)
      
      val1D <- median(as.numeric(resultsD[i,1:2]), na.rm = TRUE)
      val2D <- median(as.numeric(resultsD[i,3:4]), na.rm = TRUE)
      
      
    } else if(i %% 3 == 2){
      
      val1A <- median(as.numeric(resultsA[i,1:3]), na.rm = TRUE)
      val2A <- median(as.numeric(resultsA[i,4:6]), na.rm = TRUE)
      
      val1B <- median(as.numeric(resultsB[i,1:3]), na.rm = TRUE)
      val2B <- median(as.numeric(resultsB[i,4:6]), na.rm = TRUE)
      
      val1C <- median(as.numeric(resultsC[i,1:3]), na.rm = TRUE)
      val2C <- median(as.numeric(resultsC[i,4:6]), na.rm = TRUE)
      
      val1D <- median(as.numeric(resultsD[i,1:3]), na.rm = TRUE)
      val2D <- median(as.numeric(resultsD[i,4:6]), na.rm = TRUE)
      
    } else if(i %% 3 == 0){
      
      val1A <- median(as.numeric(resultsA[i,1:4]), na.rm = TRUE)
      val2A <- median(as.numeric(resultsA[i,5:8]), na.rm = TRUE)
      
      val1B <- median(as.numeric(resultsB[i,1:4]), na.rm = TRUE)
      val2B <- median(as.numeric(resultsB[i,5:8]), na.rm = TRUE)
      
      val1C <- median(as.numeric(resultsC[i,1:4]), na.rm = TRUE)
      val2C <- median(as.numeric(resultsC[i,5:8]), na.rm = TRUE)
      
      val1D <- median(as.numeric(resultsD[i,1:4]), na.rm = TRUE)
      val2D <- median(as.numeric(resultsD[i,5:8]), na.rm = TRUE)
      
    }
    
    vals <- c(val1A, val2A, val1B, val2B, val1C, val2C, val1D, val2D)
      
   
    str_temp <- ""
    
    if((i+1) %% 3 == 0){
      
      str_temp <- c("Gaussian ", "$\\hat{\\eta}_3$ ", "Multitaper ", "$(2\\log_2(N))^{-1}$ ",
                    "$(4\\log_2(N))^{-1}$ ", "$(8\\log_2(N))^{-1}$ ")[(i+1) %/% 3]
      
    }
    
    cat(paste(c(str_temp, c(128,512,2048)[((i-1) %% 3) + 1], vals), collapse = " & "))
    
    if(i < nrow(resultsA)){
      
      cat(" \\\\")
      
    }
    
    
    if((i %% 3) == 0 && i != 18){
      
      cat("\\hline")
      
    }
    
    cat("\n")
     
    
  }
  
}

latex_fun2(resultsA, resultsB, resultsC, resultsD)
