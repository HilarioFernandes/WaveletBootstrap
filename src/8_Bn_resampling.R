# =============================================================================
# 8_Bn_resampling.R
# =============================================================================
# Purpose  : Simulation study of Bn distribution via resampling methods.
# Chapter  : Chapter 4, Appendix C.1
# Inputs   : 7_Quasi_U_statistics_functions.R
# Outputs  : Size and power study results for the Bn statistic.
# Depends  : 7_Quasi_U_statistics_functions.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running
WORKSPACE_DIR <- file.path(BASE_PATH, "src", "WorkspaceData")
if(!dir.exists(WORKSPACE_DIR)) dir.create(WORKSPACE_DIR, recursive=TRUE)

source(file.path(BASE_PATH, "src", "7_Quasi_U_statistics_functions.R"))

# --- Testing Mode ---
TEST_MODE <- TRUE
# --------------------

# Define kernel used across all simulations
kernel <- function(x, y) (x - y)^2

n <- 100
iterations <- if(TEST_MODE) 2 else 100
B <- if(TEST_MODE) 5 else 100


################################################################################

X_e <- vector(mode = "list", length = iterations)
Y_e <- vector(mode = "list", length = iterations)

statistics_e <- NULL

statistics_boot_discy_e <- matrix(NA, nrow = iterations, ncol = B)

statistics_boot_discn_e <- matrix(NA, nrow = iterations, ncol = B)

statistics_jackknife_e <- matrix(NA, nrow = iterations, ncol = 2*n)

  
set.seed(46)

for(i in 1:iterations){
  
  X_e[[i]] <- rnorm(n)
  Y_e[[i]] <- rnorm(n)
  
  statistics_e <- c(statistics_e, B_n(X_e[[i]], Y_e[[i]], kernel))
  
}

set.seed(47)

for(i in 1:iterations){
  
  X <- X_e[[i]]
  Y <- Y_e[[i]]
  
  for(b in 1:B){
    
    #Bootstrap with discrimination
    
    temp <- bootstrap_Bn(list(X, Y), TRUE)
    
    X_boot_discy <- temp[[1]]
    Y_boot_discy <- temp[[2]]
    
    statistics_boot_discy_e[i,b] <- B_n(X_boot_discy, Y_boot_discy, kernel)
    
    #Bootstrap without discrimination
    
    temp <- bootstrap_Bn(list(X, Y), FALSE)
    
    X_boot_discn <- temp[[1]]
    Y_boot_discn <- temp[[2]]
    
    statistics_boot_discn_e[i,b] <- B_n(X_boot_discn, Y_boot_discn, kernel)
    
  }
  
}


set.seed(48)
  
for(i in 1:iterations){
  
  X <- X_e[[i]]
  Y <- Y_e[[i]]
  
  #jackknife
  
  for(m in 1:(2*n)){
    
    if(m <= n){
      
      statistics_jackknife_e[i,m] <- B_n(X[-m], Y, kernel)
      
    } else{
      
      statistics_jackknife_e[i,m] <- B_n(X, Y[-(m-n)], kernel)
      
    }
    
  }
  
}

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_e.RData"))

################################################################################

X_f <- vector(mode = "list", length = iterations)
Y_f <- vector(mode = "list", length = iterations)

statistics_f <- NULL

statistics_boot_discy_f <- matrix(NA, nrow = iterations, ncol = B)

statistics_boot_discn_f <- matrix(NA, nrow = iterations, ncol = B)

statistics_jackknife_f <- matrix(NA, nrow = iterations, ncol = 2*n)



set.seed(49)

for(i in 1:iterations){
  
  X_f[[i]] <- rnorm(n)
  Y_f[[i]] <- rnorm(n, 1, 1)
  
  statistics_f <- c(statistics_f, B_n(X_f[[i]], Y_f[[i]], kernel))
  
}



set.seed(50)

for(i in 1:iterations){
  
  X <- X_f[[i]]
  Y <- Y_f[[i]]
  
  for(b in 1:B){
    
    #Bootstrap with discrimination
    
    temp <- bootstrap_Bn(list(X, Y), TRUE)
    
    X_boot_discy <- temp[[1]]
    Y_boot_discy <- temp[[2]]
    
    statistics_boot_discy_f[i,b] <- B_n(X_boot_discy, Y_boot_discy, kernel)
    
    #Bootstrap without discrimination
    
    temp <- bootstrap_Bn(list(X, Y), FALSE)
    
    X_boot_discn <- temp[[1]]
    Y_boot_discn <- temp[[2]]
    
    statistics_boot_discn_f[i,b] <- B_n(X_boot_discn, Y_boot_discn, kernel)
    
  }
  
}



set.seed(51)

for(i in 1:iterations){
  
  X <- X_f[[i]]
  Y <- Y_f[[i]]
  
  #jackknife
  
  for(m in 1:(2*n)){
    
    if(m <= n){
      
      statistics_jackknife_f[i,m] <- B_n(X[-m], Y, kernel)
      
    } else{
      
      statistics_jackknife_f[i,m] <- B_n(X, Y[-(m-n)], kernel)
      
    }
    
  }
  
}

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_f.RData"))

################################################################################

X_g <- vector(mode = "list", length = iterations)
Y_g <- vector(mode = "list", length = iterations)

statistics_g <- NULL

statistics_boot_discy_g <- matrix(NA, nrow = iterations, ncol = B)

statistics_boot_discn_g <- matrix(NA, nrow = iterations, ncol = B)

statistics_jackknife_g <- matrix(NA, nrow = iterations, ncol = 2*n)



set.seed(52)

for(i in 1:iterations){
  
  X_g[[i]] <- rchisq(n, df = 3)
  Y_g[[i]] <- rchisq(n, df = 3)
  
  statistics_g <- c(statistics_g, B_n(X_g[[i]], Y_g[[i]], kernel))
  
}



set.seed(53)

for(i in 1:iterations){
  
  X <- X_g[[i]]
  Y <- Y_g[[i]]
  
  for(b in 1:B){
    
    #Bootstrap with discrimination
    
    temp <- bootstrap_Bn(list(X, Y), TRUE)
    
    X_boot_discy <- temp[[1]]
    Y_boot_discy <- temp[[2]]
    
    statistics_boot_discy_g[i,b] <- B_n(X_boot_discy, Y_boot_discy, kernel)
    
    #Bootstrap without discrimination
    
    temp <- bootstrap_Bn(list(X, Y), FALSE)
    
    X_boot_discn <- temp[[1]]
    Y_boot_discn <- temp[[2]]
    
    statistics_boot_discn_g[i,b] <- B_n(X_boot_discn, Y_boot_discn, kernel)
    
  }
  
}

set.seed(54)

for(i in 1:iterations){
  
  X <- X_g[[i]]
  Y <- Y_g[[i]]
  
  #jackknife
  
  for(m in 1:(2*n)){
    
    if(m <= n){
      
      statistics_jackknife_g[i,m] <- B_n(X[-m], Y, kernel)
      
    } else{
      
      statistics_jackknife_g[i,m] <- B_n(X, Y[-(m-n)], kernel)
      
    }
    
  }
  
}

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_g.RData"))

################################################################################



X_h <- vector(mode = "list", length = iterations)
Y_h <- vector(mode = "list", length = iterations)

statistics_h <- NULL

statistics_boot_discy_h <- matrix(NA, nrow = iterations, ncol = B)

statistics_boot_discn_h <- matrix(NA, nrow = iterations, ncol = B)

statistics_jackknife_h <- matrix(NA, nrow = iterations, ncol = 2*n)



set.seed(55)

for(i in 1:iterations){
  
  X_h[[i]] <- rchisq(n, df = 3)
  Y_h[[i]] <- rchisq(n, df = 6)
  
  statistics_h <- c(statistics_h, B_n(X_h[[i]], Y_h[[i]], kernel))
  
}



set.seed(56)

for(i in 1:iterations){
  
  X <- X_h[[i]]
  Y <- Y_h[[i]]
  
  for(b in 1:B){
    
    #Bootstrap with discrimination
    
    temp <- bootstrap_Bn(list(X, Y), TRUE)
    
    X_boot_discy <- temp[[1]]
    Y_boot_discy <- temp[[2]]
    
    statistics_boot_discy_h[i,b] <- B_n(X_boot_discy, Y_boot_discy, kernel)
    
    #Bootstrap without discrimination
    
    temp <- bootstrap_Bn(list(X, Y), FALSE)
    
    X_boot_discn <- temp[[1]]
    Y_boot_discn <- temp[[2]]
    
    statistics_boot_discn_h[i,b] <- B_n(X_boot_discn, Y_boot_discn, kernel)
    
  }
  
}

set.seed(57)

for(i in 1:iterations){
  
  X <- X_h[[i]]
  Y <- Y_h[[i]]
  
  #jackknife
  
  for(m in 1:(2*n)){
    
    if(m <= n){
      
      statistics_jackknife_h[i,m] <- B_n(X[-m], Y, kernel)
      
    } else{
      
      statistics_jackknife_h[i,m] <- B_n(X, Y[-(m-n)], kernel)
      
    }
    
  }
  
}


save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_h.RData"))

################################################################################

quantile_CI_boot_discy_e <- quantile_CI(statistics_boot_discy_e, 0.05)
quantile_CI_boot_discn_e <- quantile_CI(statistics_boot_discn_e, 0.05)

normal_approx_CI_boot_discy_e <- approx_CI(statistics_boot_discy_e,
                                                  "boot", 0.05)

normal_approx_CI_boot_discn_e <- approx_CI(statistics_boot_discn_e,
                                                  "boot", 0.05)

normal_approx_CI_jackknife_e <- approx_CI(statistics_jackknife_e,
                                                  "jackknife", 0.05)


quantile_CI_boot_discy_f <- quantile_CI(statistics_boot_discy_f, 0.05)
quantile_CI_boot_discn_f <- quantile_CI(statistics_boot_discn_f, 0.05)

normal_approx_CI_boot_discy_f <- approx_CI(statistics_boot_discy_f,
                                           "boot", 0.05)

normal_approx_CI_boot_discn_f <- approx_CI(statistics_boot_discn_f,
                                           "boot", 0.05)

normal_approx_CI_jackknife_f <- approx_CI(statistics_jackknife_f,
                                          "jackknife", 0.05)


quantile_CI_boot_discy_g <- quantile_CI(statistics_boot_discy_g, 0.05)
quantile_CI_boot_discn_g <- quantile_CI(statistics_boot_discn_g, 0.05)

normal_approx_CI_boot_discy_g <- approx_CI(statistics_boot_discy_g,
                                           "boot", 0.05)

normal_approx_CI_boot_discn_g <- approx_CI(statistics_boot_discn_g,
                                           "boot", 0.05)

normal_approx_CI_jackknife_g <- approx_CI(statistics_jackknife_g,
                                          "jackknife", 0.05)


quantile_CI_boot_discy_h <- quantile_CI(statistics_boot_discy_h, 0.05)
quantile_CI_boot_discn_h <- quantile_CI(statistics_boot_discn_h, 0.05)

normal_approx_CI_boot_discy_h <- approx_CI(statistics_boot_discy_h,
                                           "boot", 0.05)

normal_approx_CI_boot_discn_h <- approx_CI(statistics_boot_discn_h,
                                           "boot", 0.05)

normal_approx_CI_jackknife_h <- approx_CI(statistics_jackknife_h,
                                          "jackknife", 0.05)


################################################################################


table_results <- matrix(NA, nrow = 4, ncol = 5)
colnames(table_results) <- c("quantile_boot_disc",
                             "quantile_boot",
                             "normal_boot_disc",
                             "normal_boot",
                             "normal_jack")
rownames(table_results) <- c("e", "f", "g", "h")


table_results[1,] <- c(coverage_CI(statistics_e, t(quantile_CI_boot_discy_e)),
                       coverage_CI(statistics_e, t(quantile_CI_boot_discn_e)),
                       
                       coverage_CI(statistics_e, t(normal_approx_CI_boot_discy_e)),
                       coverage_CI(statistics_e, t(normal_approx_CI_boot_discn_e)),
                       coverage_CI(statistics_e, t(normal_approx_CI_jackknife_e)))

table_results[2,] <- c(coverage_CI(statistics_f, t(quantile_CI_boot_discy_f)),
                       coverage_CI(statistics_f, t(quantile_CI_boot_discn_f)),
                       
                       coverage_CI(statistics_f, t(normal_approx_CI_boot_discy_f)),
                       coverage_CI(statistics_f, t(normal_approx_CI_boot_discn_f)),
                       coverage_CI(statistics_f, t(normal_approx_CI_jackknife_f)))

table_results[3,] <- c(coverage_CI(statistics_g, t(quantile_CI_boot_discy_g)),
                       coverage_CI(statistics_g, t(quantile_CI_boot_discn_g)),
                       
                       coverage_CI(statistics_g, t(normal_approx_CI_boot_discy_g)),
                       coverage_CI(statistics_g, t(normal_approx_CI_boot_discn_g)),
                       coverage_CI(statistics_g, t(normal_approx_CI_jackknife_g)))

table_results[4,] <- c(coverage_CI(statistics_h, t(quantile_CI_boot_discy_h)),
                       coverage_CI(statistics_h, t(quantile_CI_boot_discn_h)),
                       
                       coverage_CI(statistics_h, t(normal_approx_CI_boot_discy_h)),
                       coverage_CI(statistics_h, t(normal_approx_CI_boot_discn_h)),
                       coverage_CI(statistics_h, t(normal_approx_CI_jackknife_h)))

View(table_results)

save.image(file.path(WORKSPACE_DIR, "8_Bn_resampling_final.RData"))

