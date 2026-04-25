# =============================================================================
# 3_Bootstrap_Wavelet_order.R
# =============================================================================
# Purpose  : Compare bw vs wb ordering and NBB vs SB bootstrap methods.
# Chapter  : Chapter 2
# Inputs   : None
# Outputs  : Distribution/variance estimation comparison outputs.
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running
WORKSPACE_DIR <- file.path(BASE_PATH, "src", "WorkspaceData")
if(!dir.exists(WORKSPACE_DIR)) dir.create(WORKSPACE_DIR, recursive=TRUE)

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

# --- Testing Mode ---
TEST_MODE <- TRUE
# --------------------

################################################################################

#Study for distribution/variance estimation

set.seed(42)

B <- if(TEST_MODE) 5 else 100
iterations <- if(TEST_MODE) 2 else 100

A.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
B.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
C.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
D.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

{

A.wv_wb_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

A.wv_bw_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

A.wv_wb_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

A.wv_bw_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

B.wv_wb_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

B.wv_bw_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

B.wv_wb_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

B.wv_bw_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

C.wv_wb_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

C.wv_bw_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

C.wv_wb_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

C.wv_bw_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

D.wv_wb_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

D.wv_bw_NBB <- list("128" = vector("list", length = iterations),
                    "512" = vector("list", length = iterations),
                    "2048" = vector("list", length = iterations))

D.wv_wb_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

D.wv_bw_SB <- list("128" = vector("list", length = iterations),
                   "512" = vector("list", length = iterations),
                   "2048" = vector("list", length = iterations))

}



# Design of the simulation study:
# 100 iterations × 100 bootstrap reps × 3 sample sizes × 4 models.

for(iter in 1:iterations){
  
  print(iter)
  
  for(i in 1:3){
    
    #print(c(iter,i))
    
    N <- 2^((2*i-1)+6)
    
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
    
    A.wv_wb_NBB[[i]][[iter]] <- bootstrap_wavelet(YA, "wb", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    A.wv_bw_NBB[[i]][[iter]] <- bootstrap_wavelet(YA, "bw", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    A.wv_wb_SB[[i]][[iter]] <- bootstrap_wavelet(YA, "wb", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    A.wv_bw_SB[[i]][[iter]] <- bootstrap_wavelet(YA, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    B.wv_wb_NBB[[i]][[iter]] <- bootstrap_wavelet(YB, "wb", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    B.wv_bw_NBB[[i]][[iter]] <- bootstrap_wavelet(YB, "bw", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    B.wv_wb_SB[[i]][[iter]] <- bootstrap_wavelet(YB, "wb", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    B.wv_bw_SB[[i]][[iter]] <- bootstrap_wavelet(YB, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    C.wv_wb_NBB[[i]][[iter]] <- bootstrap_wavelet(YC, "wb", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    C.wv_bw_NBB[[i]][[iter]] <- bootstrap_wavelet(YC, "bw", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    C.wv_wb_SB[[i]][[iter]] <- bootstrap_wavelet(YC, "wb", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    C.wv_bw_SB[[i]][[iter]] <- bootstrap_wavelet(YC, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    D.wv_wb_NBB[[i]][[iter]] <- bootstrap_wavelet(YD, "wb", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    D.wv_bw_NBB[[i]][[iter]] <- bootstrap_wavelet(YD, "bw", TRUE, "NBB", function(N){floor(4*log2(N))}, B)
    
    D.wv_wb_SB[[i]][[iter]] <- bootstrap_wavelet(YD, "wb", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
    D.wv_bw_SB[[i]][[iter]] <- bootstrap_wavelet(YD, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
    
  }
  
}

################################################################################

#Distribution comparisons

#' Calculate maximum absolute difference between empirical cdfs
#' 
#' @param list1 First nested list of wavelet variance estimates
#' @param list2 Second nested list of bootstrap wavelet variance estimates
#' @return A list of matrices (one per sample size) containing the max abs diff for each iteration and level
max_abs_diff_cdfs <- function(list1, list2){
  
  #We create a matrix in which rows correspond to iterations and columns correspond to scales.
  #Each element of this matrix is the maximum absolute difference between the bootstrap 
  #empirical cdfs
  
  list_temp <- list(matrix(NA, nrow = iterations, ncol = floor(log2(1 + (128-1)/(8-1)))),
                    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (512-1)/(8-1)))),
                    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (2048-1)/(8-1)))))
  
  
  for(i in 1:3){
    
    for(iter in 1:iterations){
      
      #we check the maximum level for which there is at least one bootstrap estimate on that particular iteration
      n_levels <- max(which(apply(list2[[i]][[iter]],2,function(x){sum(1 - is.na(x))}) >= 2))
      
      for(j in 1:n_levels){
        
        #We will then fill the matrices
        list_temp[[i]][iter, j] <- distance_cdfs(list2[[i]][[iter]][,j],list1[[i]][,j],100)
        
      }
      
    }
    
  }
  
  return(list_temp)
  
}

#comparing true values with (bootstrap then wavelet transformation)

#A.wv vs A.wv_bw_NBB
A.NBB_bw_dist <- max_abs_diff_cdfs(A.wv, A.wv_bw_NBB)

#A.wv vs A.wv_bw_SB
A.SB_bw_dist <- max_abs_diff_cdfs(A.wv, A.wv_bw_SB)

#B.wv vs B.wv_bw_NBB
B.NBB_bw_dist <- max_abs_diff_cdfs(B.wv, B.wv_bw_NBB)

#B.wv vs B.wv_bw_SB
B.SB_bw_dist <- max_abs_diff_cdfs(B.wv, B.wv_bw_SB)

#A.wv vs C.wv_bw_NBB
C.NBB_bw_dist <- max_abs_diff_cdfs(C.wv, C.wv_bw_NBB)

#C.wv vs C.wv_bw_SB
C.SB_bw_dist <- max_abs_diff_cdfs(C.wv, C.wv_bw_SB)

#D.wv vs D.wv_bw_NBB
D.NBB_bw_dist <- max_abs_diff_cdfs(D.wv, D.wv_bw_NBB)

#D.wv vs D.wv_bw_SB
D.SB_bw_dist <- max_abs_diff_cdfs(D.wv, D.wv_bw_SB)


#comparing true values with (wavelet transformation then bootstrap)

#A.wv vs A.wv_wb_NBB
A.NBB_wb_dist <- max_abs_diff_cdfs(A.wv, A.wv_wb_NBB)

#A.wv vs A.wv_wb_SB
A.SB_wb_dist <- max_abs_diff_cdfs(A.wv, A.wv_wb_SB)

#B.wv vs B.wv_wb_NBB
B.NBB_wb_dist <- max_abs_diff_cdfs(B.wv, B.wv_wb_NBB)

#B.wv vs B.wv_wb_SB
B.SB_wb_dist <- max_abs_diff_cdfs(B.wv, B.wv_wb_SB)

#A.wv vs C.wv_wb_NBB
C.NBB_wb_dist <- max_abs_diff_cdfs(C.wv, C.wv_wb_NBB)

#C.wv vs C.wv_wb_SB
C.SB_wb_dist <- max_abs_diff_cdfs(C.wv, C.wv_wb_SB)

#D.wv vs D.wv_wb_NBB
D.NBB_wb_dist <- max_abs_diff_cdfs(D.wv, D.wv_wb_NBB)

#D.wv vs D.wv_wb_SB
D.SB_wb_dist <- max_abs_diff_cdfs(D.wv, D.wv_wb_SB)


################################################################################

#Variance comparisons

#' Calculate ratio of variances
#' 
#' @param list1 First nested list of wavelet variance estimates
#' @param list2 Second nested list of bootstrap wavelet variance estimates
#' @return A list of matrices with minimum of var ratio and its inverse
ratio_vars <- function(list1, list2){
  
  #We create a matrix in which rows correspond to iterations and columns correspond to scales.
  #Each element of this matrix is the minimum between the ratio of the corresponding bootstrap
  #variance estimates and the inverse ratio
  
  list_temp <- list(matrix(NA, nrow = iterations, ncol = floor(log2(1 + (128-1)/(8-1)))),
                    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (512-1)/(8-1)))),
                    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (2048-1)/(8-1)))))
  
  
  for(i in 1:3){
    
    for(iter in 1:iterations){
      
      #we check the maximum level for which there is at least one bootstrap estimate on that particular iteration
      n_levels <- max(which(apply(list2[[i]][[iter]],2,function(x){sum(1 - is.na(x))}) >= 2))
      
      for(j in 1:n_levels){
        
        ratio <- var(list2[[i]][[iter]][,j])/var(list1[[i]][,j])
        
        #We will then fill the matrices
        list_temp[[i]][iter, j] <- min(ratio, 1/ratio)
        
      }
      
    }
    
  }
  
  return(list_temp)
  
}

#comparing true values with (bootstrap then wavelet transformation)

#A.wv vs A.wv_bw_NBB
A.NBB_bw_var <- ratio_vars(A.wv, A.wv_bw_NBB)

#A.wv vs A.wv_bw_SB
A.SB_bw_var <- ratio_vars(A.wv, A.wv_bw_SB)

#B.wv vs B.wv_bw_NBB
B.NBB_bw_var <- ratio_vars(B.wv, B.wv_bw_NBB)

#B.wv vs B.wv_bw_SB
B.SB_bw_var <- ratio_vars(B.wv, B.wv_bw_SB)

#A.wv vs C.wv_bw_NBB
C.NBB_bw_var <- ratio_vars(C.wv, C.wv_bw_NBB)

#C.wv vs C.wv_bw_SB
C.SB_bw_var <- ratio_vars(C.wv, C.wv_bw_SB)

#D.wv vs D.wv_bw_NBB
D.NBB_bw_var <- ratio_vars(D.wv, D.wv_bw_NBB)

#D.wv vs D.wv_bw_SB
D.SB_bw_var <- ratio_vars(D.wv, D.wv_bw_SB)


#comparing true values with (wavelet transformation then bootstrap)

#A.wv vs A.wv_wb_NBB
A.NBB_wb_var <- ratio_vars(A.wv, A.wv_wb_NBB)

#A.wv vs A.wv_wb_SB
A.SB_wb_var <- ratio_vars(A.wv, A.wv_wb_SB)

#B.wv vs B.wv_wb_NBB
B.NBB_wb_var <- ratio_vars(B.wv, B.wv_wb_NBB)

#B.wv vs B.wv_wb_SB
B.SB_wb_var <- ratio_vars(B.wv, B.wv_wb_SB)

#A.wv vs C.wv_wb_NBB
C.NBB_wb_var <- ratio_vars(C.wv, C.wv_wb_NBB)

#C.wv vs C.wv_wb_SB
C.SB_wb_var <- ratio_vars(C.wv, C.wv_wb_SB)

#D.wv vs D.wv_wb_NBB
D.NBB_wb_var <- ratio_vars(D.wv, D.wv_wb_NBB)

#D.wv vs D.wv_wb_SB
D.SB_wb_var <- ratio_vars(D.wv, D.wv_wb_SB)

################################################################################

#Summarizing the results (distribution estimation)

#' Summarize distribution study results
#' 
#' @param list1 A list of matrices generated by max_abs_diff_cdfs
#' @return A matrix of summarized quantiles
summary_study_dist <- function(list1){
  
  temp <- matrix(NA, nrow = 15, ncol = 8)
  
  for(i in 1:3){
    
    n_levels <- max(which(apply(list1[[i]],2,function(x){sum(1 - is.na(x))}) >= 2))
    
    for(j in 1:n_levels){
      
      temp[(5*(i-1)+1):(5*i),j] <- quantile(list1[[i]][,j], na.rm = TRUE)
      
    }
    
    
  }
  
  return(round(temp, digits = 2))
}

#bootstrap then wavelet transform

A.NBB_bw_dist_summary <- summary_study_dist(A.NBB_bw_dist)

A.SB_bw_dist_summary <- summary_study_dist(A.SB_bw_dist)

B.NBB_bw_dist_summary <- summary_study_dist(B.NBB_bw_dist)

B.SB_bw_dist_summary <- summary_study_dist(B.SB_bw_dist)

C.NBB_bw_dist_summary <- summary_study_dist(C.NBB_bw_dist)

C.SB_bw_dist_summary <- summary_study_dist(C.SB_bw_dist)

D.NBB_bw_dist_summary <- summary_study_dist(D.NBB_bw_dist)

D.SB_bw_dist_summary <- summary_study_dist(D.SB_bw_dist)

#wavelet transform then bootstrap

A.NBB_wb_dist_summary <- summary_study_dist(A.NBB_wb_dist)

A.SB_wb_dist_summary <- summary_study_dist(A.SB_wb_dist)

B.NBB_wb_dist_summary <- summary_study_dist(B.NBB_wb_dist)

B.SB_wb_dist_summary <- summary_study_dist(B.SB_wb_dist)

C.NBB_wb_dist_summary <- summary_study_dist(C.NBB_wb_dist)

C.SB_wb_dist_summary <- summary_study_dist(C.SB_wb_dist)

D.NBB_wb_dist_summary <- summary_study_dist(D.NBB_wb_dist)

D.SB_wb_dist_summary <- summary_study_dist(D.SB_wb_dist)


################################################################################


#Summarizing the results (variance estimation)

#' Summarize variance study results
#' 
#' @param list1 A list of matrices generated by ratio_vars
#' @return A matrix of summarized quantiles
summary_study_var <- function(list1){
  
  temp <- matrix(NA, nrow = 15, ncol = 8)
  
  for(i in 1:3){
    
    n_levels <- max(which(apply(list1[[i]],2,function(x){sum(1 - is.na(x))}) >= 2))
    
    for(j in 1:n_levels){
      
      temp[(5*(i-1)+1):(5*i),j] <- quantile(list1[[i]][,j], na.rm = TRUE)
      
    }
    
    
  }
  
  return(round(temp, digits = 2))
}


#bootstrap then wavelet transform

A.NBB_bw_var_summary <- summary_study_var(A.NBB_bw_var)

A.SB_bw_var_summary <- summary_study_var(A.SB_bw_var)

B.NBB_bw_var_summary <- summary_study_var(B.NBB_bw_var)

B.SB_bw_var_summary <- summary_study_var(B.SB_bw_var)

C.NBB_bw_var_summary <- summary_study_var(C.NBB_bw_var)

C.SB_bw_var_summary <- summary_study_var(C.SB_bw_var)

D.NBB_bw_var_summary <- summary_study_var(D.NBB_bw_var)

D.SB_bw_var_summary <- summary_study_var(D.SB_bw_var)

#wavelet transform then bootstrap

A.NBB_wb_var_summary <- summary_study_var(A.NBB_wb_var)

A.SB_wb_var_summary <- summary_study_var(A.SB_wb_var)

B.NBB_wb_var_summary <- summary_study_var(B.NBB_wb_var)

B.SB_wb_var_summary <- summary_study_var(B.SB_wb_var)

C.NBB_wb_var_summary <- summary_study_var(C.NBB_wb_var)

C.SB_wb_var_summary <- summary_study_var(C.SB_wb_var)

D.NBB_wb_var_summary <- summary_study_var(D.NBB_wb_var)

D.SB_wb_var_summary <- summary_study_var(D.SB_wb_var)

################################################################################

#' Generate LaTeX output for summary matrices
#' 
#' @param matrix1 First summary matrix
#' @param matrix2 Second summary matrix
#' @return Prints formatted LaTeX table rows
latex_output <- function(matrix1, matrix2){
  
  for(i in 1:15){
    
    temp1 <- ""
    temp2 <- ""
    
    if((i-1) %% 5  == 2){
      
      temp1 <- 2^((2*(1+floor((i-1)/5))-1)+6)
      
    }
    
    temp1 <- paste0(temp1, " & ", ((i-1) %% 5)*25, "\\% ")
    temp2 <- paste0(temp2, " & ")
    
    for(j in 1:(2*(2+floor((i-1)/5)))){
      
      temp1 <- paste0(temp1, " & ", matrix1[i,j])
      
      temp2 <- paste0(temp2, "& (", matrix2[i,j], ") ")
      
    }
  
    temp1 <- paste0(temp1, " \\\\ ")
    temp2 <- paste0(temp2, " \\\\ ")
    
    
    
    if((i-1) %% 5 == 4){
      
      temp2 <-paste0(temp2, " \\hline")
      
    }
    
    cat(paste0(temp1, " \n"))
    cat(paste0(temp2, " \n"))
      
  }
  
}

latex_output(A.NBB_bw_dist_summary, A.NBB_wb_dist_summary)
latex_output(A.SB_bw_dist_summary, A.SB_wb_dist_summary)

latex_output(B.NBB_bw_dist_summary, B.NBB_wb_dist_summary)
latex_output(B.SB_bw_dist_summary, B.SB_wb_dist_summary)

latex_output(C.NBB_bw_dist_summary, C.NBB_wb_dist_summary)
latex_output(C.SB_bw_dist_summary, C.SB_wb_dist_summary)

latex_output(D.NBB_bw_dist_summary, D.NBB_wb_dist_summary)
latex_output(D.SB_bw_dist_summary, D.SB_wb_dist_summary)




latex_output(A.NBB_bw_var_summary, A.NBB_wb_var_summary)
latex_output(A.SB_bw_var_summary, A.SB_wb_var_summary)

latex_output(B.NBB_bw_var_summary, B.NBB_wb_var_summary)
latex_output(B.SB_bw_var_summary, B.SB_wb_var_summary)

latex_output(C.NBB_bw_var_summary, C.NBB_wb_var_summary)
latex_output(C.SB_bw_var_summary, C.SB_wb_var_summary)

latex_output(D.NBB_bw_var_summary, D.NBB_wb_var_summary)
latex_output(D.SB_bw_var_summary, D.SB_wb_var_summary)

save.image(file.path(WORKSPACE_DIR, "3_Bootstrap_Wavelet_order.RData"))
