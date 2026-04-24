# =============================================================================
# 4_std_error.R
# =============================================================================
# Purpose  : Compare SB bootstrap SEs at 3 bandwidth choices vs multitaper estimator.
# Chapter  : Thesis Chapter 3
# Inputs   : None
# Outputs  : Standard error estimation plots (PNGs).
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

BASE_PATH <- "C:/Users/Hilar/Projects/WaveletBootstrap"  # <- SET THIS before running

source(file.path(BASE_PATH, "src", "1_Simulation_functions.R"))
source(file.path(BASE_PATH, "src", "2_Bootstrap_methods.R"))

# Set and create output directory for plots
OUTPUT_PATH <- file.path(BASE_PATH, "Plots/Plots_4")
if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)


################################################################################

#multitaper estimation

#install.packages("multitaper")
library(multitaper)

#' Calculate the wavelet variance standard deviation estimates using multitaper
#' 
#' @param Z List with all squared wavelet coefficients (each element corresponds to a scale)
#' @param v List with all tapers (each element corresponds to a scale)
#' @return A vector of multitaper standard error estimates
sd_estimate_multitaper <- function(Z, v){
  
  n_levels <- length(Z)
  
  N <- nrow(v[[1]])
  
  L <- (2^(1:n_levels) -1)*(8 - 1)+1
  
  #now we calculate the computable quantity J(0) and V(0)
  
  J <- matrix(0,ncol = n_levels, nrow = 5)
  V <- matrix(0,ncol = n_levels, nrow = 5)
  
  for(k in 1:5){
    
    for(j in 1:n_levels){
      
      #using non-boundary coeffs and the corresponding tapers
      #must set remove_boundary_coeffs <- TRUE
      J[k,j] <- Z[[j]] %*% v[[j]][L[j]:N,k]
      
      
      #V[k,j] <- sum(v[[j]][L[j]:N,k])
      #V[k,j] <- sum(v[[j]][1:(N-L[j]+1),k])
      V[k,j] <- sum(v[[j]][,k]) 
      
    }
    
  }
  
  #calculating the estimated variance
  
  estimates <- NULL
  
  for(j in 1:n_levels){
    
    #denominator of (12)
    
    denominator <- V[1,j]^2 + V[3,j]^2 + V[5,j]^2
    
    temp <- (J[1,j]*V[1,j] + J[3,j]*V[3,j] + J[5,j]*V[5,j])/denominator
    
    sum <- 0
    
    for(k in 1:5){
      
      sum <- sum + ((J[k,j] - temp*V[k,j])^2)/(N-L[j]+1)
      
    }
    
    sum <- sum/5
    
    estimates <- c(estimates, sum)
    
  }
  
  return(sqrt(estimates))
  
}

#' Generate discrete prolate spheroidal sequences (DPSS) taper
#'
#' @param N Sample size
#' @param L Filter length for the current scale
#' @return DPSS taper vector
taper_fun <- function(N,L){
  
  #return(dpss(N,5,3.5/(N-L+1))$v)
  
  #return(dpss(N,5,7/(N-L+1))$v)
  
  return(dpss(N,5,7)$v)
  
}

#given a vector of filter lengths L <- (2^(1:n_levels) -1)*(8 - 1)+1,
#calculate v <- lapply(L,taper_fun) and use it on sd_estimate_multitaper()

################################################################################

#study for standard error
# Simulation design: 4 models × 3 sample sizes × 3 bandwidths

set.seed(43)

B <- 100
iterations <- 100

A.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
B.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
C.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
D.wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

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

{
  
  A.wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  B.wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  C.wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  D.wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
  
}


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
    
    #multitaper wavelet variance standard errors (via multitaper)
    
    n_levels <- floor(log2(1+(N-1)/(8-1)))
    
    L <- (2^(1:n_levels) -1)*(8 - 1)+1
    
    taper_fun_temp <- function(L){
      return(taper_fun(N,L))
    }
    
    v <- lapply(L,taper_fun_temp)
    
    A.wv_sd_multitaper[[i]] <- rbind(A.wv_sd_multitaper[[i]], sd_estimate_multitaper(nonboundary_squared_MODWT_coeffs(YA),v))
    B.wv_sd_multitaper[[i]] <- rbind(B.wv_sd_multitaper[[i]], sd_estimate_multitaper(nonboundary_squared_MODWT_coeffs(YB),v))
    C.wv_sd_multitaper[[i]] <- rbind(C.wv_sd_multitaper[[i]], sd_estimate_multitaper(nonboundary_squared_MODWT_coeffs(YC),v))
    D.wv_sd_multitaper[[i]] <- rbind(D.wv_sd_multitaper[[i]], sd_estimate_multitaper(nonboundary_squared_MODWT_coeffs(YD),v))
    
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

#standard errors

std_true <- function(list1){
  
  temp <- list(NULL,NULL,NULL)
  
  for(i in 1:3){
    
    temp[[i]] <- apply(list1[[i]], 2, function(x){sd(x, na.rm = TRUE)})
    
  }
  
  return(temp)
  
}


A.wv_std <- std_true(A.wv)
B.wv_std <- std_true(B.wv)
C.wv_std <- std_true(C.wv)
D.wv_std <- std_true(D.wv)

std_boot <- function(list1){
  
  temp <- list(NULL,NULL,NULL)
  
  for(i in 1:3){
    
    for(iteration in 1:iterations){
      
      temp[[i]] <- rbind(temp[[i]], apply(list1[[i]][[iteration]], 2, function(x){sd(x)}))
      
    }
    
  }
  
  return(temp)
  
}


A.wv_SB_2_std <- std_boot(A.wv_SB_2)
A.wv_SB_4_std <- std_boot(A.wv_SB_4)
A.wv_SB_8_std <- std_boot(A.wv_SB_8)

B.wv_SB_2_std <- std_boot(B.wv_SB_2)
B.wv_SB_4_std <- std_boot(B.wv_SB_4)
B.wv_SB_8_std <- std_boot(B.wv_SB_8)

C.wv_SB_2_std <- std_boot(C.wv_SB_2)
C.wv_SB_4_std <- std_boot(C.wv_SB_4)
C.wv_SB_8_std <- std_boot(C.wv_SB_8)

D.wv_SB_2_std <- std_boot(D.wv_SB_2)
D.wv_SB_4_std <- std_boot(D.wv_SB_4)
D.wv_SB_8_std <- std_boot(D.wv_SB_8)





################################################################################

#plots

par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))

col1 <- "darkgrey"
col2 <- "lightgrey"


## ---- Model A ----

{
  png(file=file.path(OUTPUT_PATH, "A_wv_std_2.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(A.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.17), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(A.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(A.wv_std[[1]],
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()

  plot(A.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(A.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[2]],
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()

  plot(A.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(A.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[3]],
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "A_wv_std_4.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(A.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.17), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(A.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(A.wv_std[[1]],
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()

  plot(A.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(A.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[2]],
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)


  plot(A.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(A.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[3]],
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "A_wv_std_8.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(A.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.17), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(A.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(A.wv_std[[1]],
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  
  plot(A.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(A.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[2]],
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(A.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(A.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(A.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(A.wv_std[[3]],
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(A.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
}


## ---- Model B ----

{
  png(file=file.path(OUTPUT_PATH, "B_wv_std_2.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(B.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.5), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(B.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(B.wv_std[[1]],
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()
  
  plot(B.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(B.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[2]],
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()
  
  plot(B.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(B.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[3]],
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "B_wv_std_4.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(B.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.5), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(B.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(B.wv_std[[1]],
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  #dev.off()
  
  plot(B.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(B.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[2]],
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  
  plot(B.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(B.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[3]],
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "B_wv_std_8.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(B.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.5), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(B.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(B.wv_std[[1]],
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  
  plot(B.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(B.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[2]],
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  
  plot(B.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(B.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(B.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(B.wv_std[[3]],
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(B.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
}

## ---- Model C ----

{
  
  png(file=file.path(OUTPUT_PATH, "C_wv_std_2.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(C.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.18), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(C.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(C.wv_std[[1]],
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(C.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[2]],
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(C.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[3]],
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "C_wv_std_4.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(C.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", xaxt = "n", ylim = c(0,0.18),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(C.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(C.wv_std[[1]],
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(C.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[2]],
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(C.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[3]],
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "C_wv_std_8.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(C.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", xaxt = "n", ylim = c(0,0.18),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(C.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(C.wv_std[[1]],
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(C.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[2]],
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(C.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(C.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(C.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(C.wv_std[[3]],
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(C.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
}

## ---- Model D ----

{
  
  png(file=file.path(OUTPUT_PATH, "D_wv_std_2.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(D.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.12), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(D.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(D.wv_std[[1]],
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)", ylim = c(0,0.07),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(D.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[2]],
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)", ylim = c(0,0.03),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(D.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[3]],
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_2_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "D_wv_std_4.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(D.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.12), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(D.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(D.wv_std[[1]],
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)", ylim = c(0,0.07),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(D.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[2]],
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)", ylim = c(0,0.03),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(D.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[3]],
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_4_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
  #
  
  png(file=file.path(OUTPUT_PATH, "D_wv_std_8.png"), width=1800, height=2100, res=210)
  par(mfrow=c(3,1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(D.wv_std[[1]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=128)", ylim = c(0,0.12), xaxt = "n",
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  axis(1, at= 1:4, cex.axis=1.5)
  polygon(c(1:4, 4:1), c(apply(D.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  polygon(c(1:4, 4:1), c(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  lines(D.wv_std[[1]],
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[1]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[1]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[2]], type = "l", xlab = "", ylab = "Desvio padrão \n (N=512)", ylim = c(0,0.07),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5, xlim = c(1,8))
  polygon(c(1:6, 6:1), c(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:6, 6:1), c(apply(D.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[2]],
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[2]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[2]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)

  plot(D.wv_std[[3]], type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (N=2048)", ylim = c(0,0.030),
       cex.lab=1.5, cex.main=1.5, cex.axis=1.5, cex.sub=1.5,
       lwd = 1.5)
  polygon(c(1:8, 8:1), c(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col1, border = "NA")
  polygon(c(1:8, 8:1), c(apply(D.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}),
                         rev(apply(D.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}))),
          col = col2, border = "NA")
  lines(D.wv_std[[3]],
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[2]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_SB_8_std[[3]], 2, function(x){quantile(x)[[4]]}), lty = "dotted",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[2]]}), lty = "longdash",
        lwd = 1.5)
  lines(apply(D.wv_sd_multitaper[[3]], 2, function(x){quantile(x)[[4]]}), lty = "longdash",
        lwd = 1.5)
  legend("topleft", legend=c("Multitaper", "Bootstrap"),  
         fill = c(col1, col2), cex = 1.5)
  dev.off()
  
}
