# =============================================================================
# 11_Changepoint_sim.R
# =============================================================================
# Purpose  : Simulation study for changepoint detection performance.
# Chapter  : Thesis Chapter 3
# Inputs   : 10_Comparison_Changepoint_functions.R
# Outputs  : Power and size rejection rates for changepoint detection.
# Depends  : 10_Comparison_Changepoint_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "10_Comparison_Changepoint_functions.R"))

################################################################################

#' Run changepoint detection simulation (N=2048, block=128)
#' 
#' @param model Model identifier ("A" through "D")
#' @param iterations Number of simulation iterations
#' @param B Number of bootstrap reps
#' @param alpha_0 Target significance level
#' @param change Logical. If TRUE, introduces a variance change.
#' @param change_type Either "abrupt" or "drift"
#' @return Matrix of rejection proportions per block
changepoint_sim <- function(model, iterations, B, alpha_0, change, change_type){
  
  sim_fun <- get(paste0("Model_", model, "_sim"))
  
  results <- matrix(0, nrow = 15, ncol = 4)
  
  rownames(results) <- 1:15
  
  for(i in 1:iterations){
    
    print(i)
    
    X <- sim_fun(2048)
    
    if(change == TRUE){
      
      if(change_type == "abrupt"){
        
        X[501:1471] <- 1.25*X[501:1471]
        X[1472:2048] <- (1.25^2)*X[1472:2048]
        
      } else if (change_type == "drift"){
        
        X <- seq(1,1.25^2, length.out = length(X))*X
        
      }
      
    }
    
    #plot(X, type = "l")
    
    temp <- changepointdetection(X, 128, FALSE, "F", alpha_0, 2, NA)
    results[temp[,2]-1,1] <- results[temp[,2]-1,1] + 1
    
    temp <- changepointdetection(X, 128, TRUE, "F", alpha_0, 2, NA)
    results[temp[,2]-1,2] <- results[temp[,2]-1,2] + 1
    
    temp <- changepointdetection(X, 128, FALSE, "boot_quant", alpha_0, 2, B)
    results[temp[,2]-1,3] <- results[temp[,2]-1,3] + 1
    
    temp <- changepointdetection(X, 128, TRUE, "boot_quant", alpha_0, 2, B)
    results[temp[,2]-1,4] <- results[temp[,2]-1,4] + 1
  
  }
  
  return(results/iterations)
  
}

changepoint_sim_B_abrupt <- changepoint_sim("B", 10, 100, 0.05, TRUE, "abrupt")
changepoint_sim_B_drift <- changepoint_sim("B", 10, 100, 0.05, TRUE, "drift")
changepoint_sim_B_NA <- changepoint_sim("B", 10, 100, 0.05, FALSE, NA)

################################################################################

iterations <- 100
B <- 100

################################################################################
#128

set.seed(66)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_A_abrupt <- changepoint_sim("A", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_A_drift <- changepoint_sim("A", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_A_NA <- changepoint_sim("A", iterations, B, 0.05, FALSE, NA)
  
  changepoint_sim_B_abrupt <- changepoint_sim("B", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_B_drift <- changepoint_sim("B", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_B_NA <- changepoint_sim("B", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

set.seed(67)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_C_abrupt <- changepoint_sim("C", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_C_drift <- changepoint_sim("C", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_C_NA <- changepoint_sim("C", iterations, B, 0.05, FALSE, NA)
  
  changepoint_sim_D_abrupt <- changepoint_sim("D", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_D_drift <- changepoint_sim("D", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_D_NA <- changepoint_sim("D", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

################################################################################

#' Run changepoint detection simulation (N=8192, block=512)
#' 
#' @param model Model identifier ("A" through "D")
#' @param iterations Number of simulation iterations
#' @param B Number of bootstrap reps
#' @param alpha_0 Target significance level
#' @param change Logical. If TRUE, introduces a variance change.
#' @param change_type Either "abrupt" or "drift"
#' @return Matrix of rejection proportions per block
changepoint_sim_512 <- function(model, iterations, B, alpha_0, change, change_type){
  
  sim_fun <- get(paste0("Model_", model, "_sim"))
  
  results <- matrix(0, nrow = 15, ncol = 4)
  
  rownames(results) <- 1:15
  
  for(i in 1:iterations){
    
    print(i)
    
    X <- sim_fun(8192)
    
    if(change == TRUE){
      
      if(change_type == "abrupt"){
        
        X[2004:5884] <- 1.25*X[2004:5884]
        X[5885:8192] <- (1.25^2)*X[5885:8192]
        
      } else if (change_type == "drift"){
        
        X <- seq(1,1.25^2, length.out = length(X))*X
        
      }
      
    }
    
    #plot(X, type = "l")
    
    temp <- changepointdetection(X, 512, FALSE, "F", alpha_0, 3, NA)
    results[temp[,2]-1,1] <- results[temp[,2]-1,1] + 1
    
    temp <- changepointdetection(X, 512, TRUE, "F", alpha_0, 3, NA)
    results[temp[,2]-1,2] <- results[temp[,2]-1,2] + 1
    
    temp <- changepointdetection(X, 512, FALSE, "boot_quant", alpha_0, 3, B)
    results[temp[,2]-1,3] <- results[temp[,2]-1,3] + 1
    
    temp <- changepointdetection(X, 512, TRUE, "boot_quant", alpha_0, 3, B)
    results[temp[,2]-1,4] <- results[temp[,2]-1,4] + 1
    
  }
  
  return(results/iterations)
  
}

changepoint_sim_512("B", 10, 100, 0.05, TRUE, "abrupt")

################################################################################
#512

set.seed(68)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_A_abrupt_512 <- changepoint_sim_512("A", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_A_drift_512 <- changepoint_sim_512("A", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_A_NA_512 <- changepoint_sim_512("A", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

set.seed(69)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_B_abrupt_512 <- changepoint_sim_512("B", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_B_drift_512 <- changepoint_sim_512("B", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_B_NA_512 <- changepoint_sim_512("B", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

set.seed(70)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_C_abrupt_512 <- changepoint_sim_512("C", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_C_drift_512 <- changepoint_sim_512("C", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_C_NA_512 <- changepoint_sim_512("C", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

set.seed(71)

{
  
  start_time <- Sys.time()
  
  changepoint_sim_D_abrupt_512 <- changepoint_sim_512("D", iterations, B, 0.05, TRUE, "abrupt")
  changepoint_sim_D_drift_512 <- changepoint_sim_512("D", iterations, B, 0.05, TRUE, "drift")
  changepoint_sim_D_NA_512 <- changepoint_sim_512("D", iterations, B, 0.05, FALSE, NA)
  
  end_time <- Sys.time()
  
  duration <- end_time - start_time
  print(duration)
  
}

################################################################################

#' Generate LaTeX table for changepoint simulation results
#' 
#' @param df Results matrix for first sample size (e.g., 128)
#' @param df2 Results matrix for second sample size (e.g., 512)
#' @return Prints a formatted LaTeX table to the console
latex_fun <- function(df, df2){
  
  cat("\\begin{table}[h!] \n")
  
  cat("\\centering \n")
  
  cat("\\caption{Taxas empíricas de rejeição para $l=$ e modelo () com/sem mudança abrupta/gradual de variância.} \n")
  
  cat("\\label{temp} \n")
  
  cat("\\begin{tabular}{c|cccc} \n")
  
  cat("Update & Aprox. $F$ & Aprox. $F$ (seq.) & Bootstrap & Bootstrap (seq.) \\\\ \\hline \n")
  
  for(i in 1:15){
    
    cat(paste0(i, " & ", paste(df[i,], collapse = " & "), " \\\\ \n"))
    cat(paste0(" & (", paste(df2[i,], collapse = ") & ("), ") \\\\ \\hline \n"))
    
  }
  
  cat("\\end{tabular} \n")
  
  cat("\\end{table} \n")
  
}

latex_fun(changepoint_sim_A_NA, changepoint_sim_A_NA_512)
latex_fun(changepoint_sim_A_drift, changepoint_sim_A_drift_512)
latex_fun(changepoint_sim_A_abrupt, changepoint_sim_A_abrupt_512)

latex_fun(changepoint_sim_B_NA, changepoint_sim_B_NA_512)
latex_fun(changepoint_sim_B_drift, changepoint_sim_B_drift_512)
latex_fun(changepoint_sim_B_abrupt, changepoint_sim_B_abrupt_512)

latex_fun(changepoint_sim_C_NA, changepoint_sim_C_NA_512)
latex_fun(changepoint_sim_C_drift, changepoint_sim_C_drift_512)
latex_fun(changepoint_sim_C_abrupt, changepoint_sim_C_abrupt_512)

latex_fun(changepoint_sim_D_NA, changepoint_sim_D_NA_512)
latex_fun(changepoint_sim_D_drift, changepoint_sim_D_drift_512)
latex_fun(changepoint_sim_D_abrupt, changepoint_sim_D_abrupt_512)

