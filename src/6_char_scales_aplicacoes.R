# =============================================================================
# 6_char_scales_aplicacoes.R
# =============================================================================
# Purpose  : Real-data characteristic scales applications (MJO, INMET datasets).
# Chapter  : Thesis Chapter 3
# Inputs   : MJO.txt, dados_INMET_processados.csv (Data files not included).
# Outputs  : Characteristic scale estimates and application plots.
# Depends  : 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

library(imputeTS)

source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

B <- 100

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

#' beta_hat calculation for quadratic interpolation
#' 
#' @param wv_est Vector of wavelet variance estimates
#' @param j_0 Target level index
#' @return Vector of beta coefficients
beta_hat_est <- function(wv_est, j_0){
  
  H <- rbind(c(-1/2,0,1/2),c(1,-2,1)) 
  
  return(as.vector(H %*% log2(wv_est[(j_0 - 1):(j_0 + 1)])))
  
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

#' Function that calculates the var-cov matrix associated to \hat{\beta}
#' 
#' @param sigma2_est The Sigma2 matrix estimate
#' @return The variance-covariance matrix for beta coefficients
Var_beta_est <- function(sigma2_est){
  
  H <- rbind(c(-1/2,0,1/2),c(1,-2,1)) 
  
  return(H %*% sigma2_est %*% t(H))
  
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

#' Calculate bootstrap confidence intervals for characteristic scales
#' 
#' @param char_scale_est_boot Vector of bootstrap characteristic scale estimates
#' @param alpha Significance level
#' @param B Number of bootstrap reps
#' @return A vector with [lower, upper] limits
CI_calc_boot <- function(char_scale_est_boot, alpha, B){
  
  output <- list(c(NULL,NULL),c(NULL,NULL),c(NULL,NULL))

  if(length(which(char_scale_est_boot == -1)) > B/2){
        
      return(c(-1,-1))
        
  } else{
    
    x <- char_scale_est_boot
        
    return(quantile(x[x > 0], c(alpha/2,1 - alpha/2)))
        
  }

}

################################################################################

#' Wrapper for characteristic scale estimation and CI calculation
#' 
#' @param data Time series data vector
#' @param alpha Significance level
#' @param B Number of bootstrap reps
#' @return A list containing point estimate, literature CI and bootstrap CI
char_scale_fun <- function(data, alpha, B){
  
  wv_est <- wv_estimates(data)
  
  wv_boot <- NULL
  
  N <- length(data)
  
  n_levels <- floor(log2(1+(N-1)/(8-1)))
  
  wv_est <- wv_estimates(data)
  
  char_scale_est <- Char_scale_est(wv_est,max_char_scale(wv_est))
  
  wv_cs_CI <- CI_calc(data, alpha) 
  
  wv_boot <- bootstrap_wavelet(data, "bw", TRUE, "SB", function(N){1/(4*log2(N))}, B)
  
  char_scale_est_boot <- NULL
  
  for(i in 1:B){
    
    x <- wv_boot[i,]
    
    char_scale_est_boot <- c(char_scale_est_boot, Char_scale_est(x, max_char_scale(x)))
    
  }
  
  wv_cs_CI_boot <- CI_calc_boot(char_scale_est_boot, 0.05, B)
  
  names(wv_cs_CI) <- c("2.5%", "97.5%")
  
  return(list("Point_estimate" = unname(char_scale_est),
              "Literature" = wv_cs_CI,
              "Bootstrap" = wv_cs_CI_boot))
  
}

################################################################################

#MJO

# Note: Data file is not included in the repository. See thesis for sources.
MJO <- read.table(file.path(BASE_PATH, "Dados", "MJO", "MJO.txt"))

MJO <- as.numeric(MJO$INDEX_1[2:nrow(MJO)])

MJO <- na_interpolation(MJO)

plot(MJO, type = "l")

MJO_paper <- MJO[1:2354]

plot(MJO_paper, type = "l")

set.seed(0)

results_MJO <- char_scale_fun(MJO_paper, 0.05, 100)

{
  png(file=file.path(OUTPUT_PATH, "plot_aplicacao.png"), width=1800, height=900, res = 210)
  
  par(mfrow=c(1,1))
  par(mar=c(4.1, 4.1, 1, 2.1))
  
  wv_est <- wv_estimates(MJO_paper)
  
  plot(log2(wv_est), xlab = "Nível", ylab = "Var. de ondaletas (log)", pch = 19,
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5, ylim = c(-7,-1))
  
  x <- 2:4
  y <- log2(wv_est)[2:4]
  quadratic <- solve(cbind(1, x, x^2), y)
  
  quadratic_fun <- function(x){
    
    quadratic[[1]] + quadratic[[2]]*x + quadratic[[3]]*x^2
    
  }
  
  lines(seq(2,4,length.out = 50), quadratic_fun(seq(2,4,length.out = 50)),
        lwd = 1.5)
  
  abline(v = results_MJO[[1]], lty = "dashed",
         lwd = 1.5)
  
  dev.off()
}

################################################################################

#INMET data

# Note: Data file is not included in the repository. See thesis for sources.
INMET <- read.csv(file.path(BASE_PATH, "Dados", "Selecao", "dados_INMET_processados.csv"))

results_INMET <- matrix(NA, nrow = 10, ncol = 6)

set.seed(0)

for(i in 1:10){
  print(i)
  
  temp <- char_scale_fun(INMET[,i+2], 0.05, B)
  
  results_INMET[i,1] <- colnames(INMET)[i+2]
  results_INMET[i,-1] <- unlist(temp)
  
}



{
  png(file=file.path(OUTPUT_PATH, "plot_aplicacao_INMET.png"), width=1800, height=1200, res = 210)
  
par(mfrow=c(2,2))
par(mar=c(4,5,1,1))

wv_est <- wv_estimates(INMET[,3])
plot(log2(wv_est), xlab = "", ylab = "Var. de ondaletas (log)", pch = 19, xaxt = "n", ylim = c(-25,0),
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5)
Axis(side=1, labels=FALSE)

x <- 3:5
y <- log2(wv_est)[3:5]
quadratic <- solve(cbind(1, x, x^2), y)

quadratic_fun <- function(x){
  
  quadratic[[1]] + quadratic[[2]]*x + quadratic[[3]]*x^2
  
}

lines(seq(3,5,length.out = 50), quadratic_fun(seq(3,5,length.out = 50)),
      lwd = 1.5)

abline(v = results_INMET[1,2], lty = "dashed",
       lwd = 1.5)

wv_est <- wv_estimates(INMET[,4])
plot(log2(wv_est), xlab = "", ylab = "", pch = 19, xaxt = "n",ylim = c(-25,0),
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5)
Axis(side=1, labels=FALSE)

wv_est <- wv_estimates(INMET[,5])
plot(log2(wv_est), xlab = "Nível", ylab = "Var. de ondaletas (log)", pch = 19,ylim = c(-25,0),
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5)

wv_est <- wv_estimates(INMET[,6])
plot(log2(wv_est), xlab = "Nível", ylab = "", pch = 19,ylim = c(-25,0),
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5)

dev.off()
}

wv_ests <- NULL

for(i in 3:12){
  
  wv_ests <- rbind(wv_ests, wv_estimates(INMET[,i]))
  
}
