# =============================================================================
# 4_std_error.R
# =============================================================================
# Purpose  : Compare SB bootstrap SEs at 3 bandwidth choices vs multitaper estimator.
# Chapter  : Chapter 3, Appendix b.1
# Inputs   : None
# Outputs  : Standard error estimation plots (PNGs).
# Depends  : 1_Simulation_functions.R, 2_Bootstrap_methods.R
# Author   : Hilário Fernandes de Araujo Júnior
# Date     : 2024
# =============================================================================

base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap" # <- SET THIS before running
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "1_Simulation_functions.R"))
source(file.path(base_path, "src", "2_Bootstrap_methods.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

# Set and create output directory for plots
output_path <- file.path(base_path, "Plots/Plots_plots_4")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


################################################################################

# multitaper estimation

# install.packages("multitaper")
library(multitaper)

#' Calculate the wavelet variance standard deviation estimates using multitaper
#'
#' @param z List with all squared wavelet coefficients (each element corresponds to a scale)
#' @param v List with all tapers (each element corresponds to a scale)
#' @return A vector of multitaper standard error estimates
sd_estimate_multitaper <- function(z, v) {
  n_levels <- length(z)

  n <- nrow(v[[1]])

  l <- (2^(1:n_levels) - 1) * (8 - 1) + 1

  # now we calculate the computable quantity J(0) and V(0)

  J <- matrix(0, ncol = n_levels, nrow = 5)
  V <- matrix(0, ncol = n_levels, nrow = 5)

  for (k in 1:5) {
    for (j in 1:n_levels) {
      # using non-boundary coeffs and the corresponding tapers
      # must set remove_boundary_coeffs <- TRUE
      J[k, j] <- z[[j]] %*% v[[j]][l[j]:n, k]


      # V[k,j] <- sum(v[[j]][l[j]:n,k])
      # V[k,j] <- sum(v[[j]][1:(n-l[j]+1),k])
      V[k, j] <- sum(v[[j]][, k])
    }
  }

  # calculating the estimated variance

  estimates <- NULL

  for (j in 1:n_levels) {
    # denominator of (12)

    denominator <- V[1, j]^2 + V[3, j]^2 + V[5, j]^2

    temp <- (J[1, j] * V[1, j] + J[3, j] * V[3, j] + J[5, j] * V[5, j]) / denominator

    sum <- 0

    for (k in 1:5) {
      sum <- sum + ((J[k, j] - temp * V[k, j])^2) / (n - l[j] + 1)
    }

    sum <- sum / 5

    estimates <- c(estimates, sum)
  }

  return(sqrt(estimates))
}

#' Generate discrete prolate spheroidal sequences (DPSS) taper
#'
#' @param n Sample size
#' @param l Filter length for the current scale
#' @return DPSS taper vector
taper_fun <- function(n, l) {
  # return(dpss(n,5,3.5/(n-l+1))$v)

  # return(dpss(n,5,7/(n-l+1))$v)

  return(dpss(n, 5, 7)$v)
}

# given a vector of filter lengths l <- (2^(1:n_levels) -1)*(8 - 1)+1,
# calculate v <- lapply(l,taper_fun) and use it on sd_estimate_multitaper()

################################################################################

# study for standard error
# Simulation design: 4 models × 3 sample sizes × 3 bandwidths

set.seed(43)

b <- if (test_mode) 5 else 100
iterations <- if (test_mode) 2 else 100

a_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
b_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
c_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
d_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

a_wv_sb_2 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_sb_4 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_sb_8 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_sb_2 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_sb_4 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_sb_8 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_sb_2 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_sb_4 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_sb_8 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_sb_2 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_sb_4 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_sb_8 <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
b_wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
c_wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)
d_wv_sd_multitaper <- list("128" = NULL, "512" = NULL, "2048" = NULL)


for (iter in 1:iterations) {
  print(iter)

  for (i in 1:3) {
    # print(c(iter,i))

    n <- 2^((2 * i - 1) + 6)

    # simulating from each model
    ya <- model_a_sim(n)
    yb <- model_b_sim(n)
    yc <- model_c_sim(n)
    yd <- model_d_sim(n)

    # calculating the point estimates of the wavelet variances
    a_wv[[i]] <- rbind(a_wv[[i]], wv_estimates(ya))
    b_wv[[i]] <- rbind(b_wv[[i]], wv_estimates(yb))
    c_wv[[i]] <- rbind(c_wv[[i]], wv_estimates(yc))
    d_wv[[i]] <- rbind(d_wv[[i]], wv_estimates(yd))

    # multitaper wavelet variance standard errors (via multitaper)

    n_levels <- floor(log2(1 + (n - 1) / (8 - 1)))

    l <- (2^(1:n_levels) - 1) * (8 - 1) + 1

    taper_fun_temp <- function(l) {
      return(taper_fun(n, l))
    }

    v <- lapply(l, taper_fun_temp)

    a_wv_sd_multitaper[[i]] <- rbind(a_wv_sd_multitaper[[i]], sd_estimate_multitaper(nonbnd_sq_modwt_coeffs(ya), v))
    b_wv_sd_multitaper[[i]] <- rbind(b_wv_sd_multitaper[[i]], sd_estimate_multitaper(nonbnd_sq_modwt_coeffs(yb), v))
    c_wv_sd_multitaper[[i]] <- rbind(c_wv_sd_multitaper[[i]], sd_estimate_multitaper(nonbnd_sq_modwt_coeffs(yc), v))
    d_wv_sd_multitaper[[i]] <- rbind(d_wv_sd_multitaper[[i]], sd_estimate_multitaper(nonbnd_sq_modwt_coeffs(yd), v))

    # bootstrap wavelet estimates

    a_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    a_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    a_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(ya, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    b_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    b_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    b_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yb, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    c_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    c_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    c_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yc, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)

    d_wv_sb_2[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (2 * log2(n))
    }, b)

    d_wv_sb_4[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (4 * log2(n))
    }, b)

    d_wv_sb_8[[i]][[iter]] <- bootstrap_wavelet(yd, "bw", TRUE, "SB", function(n) {
      1 / (8 * log2(n))
    }, b)
  }
}

################################################################################

# standard errors

std_true <- function(list1) {
  temp <- list(NULL, NULL, NULL)

  for (i in 1:3) {
    temp[[i]] <- apply(list1[[i]], 2, function(x) {
      sd(x, na.rm = TRUE)
    })
  }

  return(temp)
}


a_wv_std <- std_true(a_wv)
b_wv_std <- std_true(b_wv)
c_wv_std <- std_true(c_wv)
d_wv_std <- std_true(d_wv)

std_boot <- function(list1) {
  temp <- list(NULL, NULL, NULL)

  for (i in 1:3) {
    for (iteration in 1:iterations) {
      temp[[i]] <- rbind(temp[[i]], apply(list1[[i]][[iteration]], 2, function(x) {
        sd(x)
      }))
    }
  }

  return(temp)
}


a_wv_sb_2_std <- std_boot(a_wv_sb_2)
a_wv_sb_4_std <- std_boot(a_wv_sb_4)
a_wv_sb_8_std <- std_boot(a_wv_sb_8)

b_wv_sb_2_std <- std_boot(b_wv_sb_2)
b_wv_sb_4_std <- std_boot(b_wv_sb_4)
b_wv_sb_8_std <- std_boot(b_wv_sb_8)

c_wv_sb_2_std <- std_boot(c_wv_sb_2)
c_wv_sb_4_std <- std_boot(c_wv_sb_4)
c_wv_sb_8_std <- std_boot(c_wv_sb_8)

d_wv_sb_2_std <- std_boot(d_wv_sb_2)
d_wv_sb_4_std <- std_boot(d_wv_sb_4)
d_wv_sb_8_std <- std_boot(d_wv_sb_8)


################################################################################

# plots

par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))

col1 <- "darkgrey"
col2 <- "lightgrey"


## ---- Model A ----

{
  png(file = file.path(output_path, "A_wv_std_2.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(a_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.17), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(a_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(a_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(a_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "A_wv_std_4.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(a_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.17), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(a_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(a_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )


  plot(a_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "A_wv_std_8.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(a_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.17), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(a_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(a_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(a_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(a_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(a_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(a_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(a_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(a_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()
}


## ---- Model b ----

{
  png(file = file.path(output_path, "B_wv_std_2.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(b_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.5), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(b_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(b_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(b_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "B_wv_std_4.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(b_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.5), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(b_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  # dev.off()

  plot(b_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(b_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "B_wv_std_8.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(b_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.5), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(b_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(b_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(b_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(b_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(b_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(b_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(b_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(b_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()
}

## ---- Model C ----

{
  png(file = file.path(output_path, "C_wv_std_2.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(c_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.18), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(c_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "C_wv_std_4.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(c_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", xaxt = "n", ylim = c(0, 0.18),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(c_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "C_wv_std_8.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(c_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", xaxt = "n", ylim = c(0, 0.18),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(c_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(c_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(c_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(c_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(c_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(c_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(c_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()
}

## ---- Model D ----

{
  png(file = file.path(output_path, "D_wv_std_2.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(d_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.12), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(d_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)", ylim = c(0, 0.07),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)", ylim = c(0, 0.03),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_2_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "D_wv_std_4.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(d_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.12), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(d_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)", ylim = c(0, 0.07),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)", ylim = c(0, 0.03),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_4_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()

  #

  png(file = file.path(output_path, "D_wv_std_8.png"), width = 1800, height = 2100, res = 210)
  par(mfrow = c(3, 1), mar = c(4.1, 6.1, 1.1, 2.1))
  plot(d_wv_std[[1]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=128)", ylim = c(0, 0.12), xaxt = "n",
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  axis(1, at = 1:4, cex.axis = 1.5)
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  polygon(c(1:4, 4:1), c(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  lines(d_wv_std[[1]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[1]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[2]],
    type = "l", xlab = "", ylab = "Desvio padrão \n (n=512)", ylim = c(0, 0.07),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5, xlim = c(1, 8)
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:6, 6:1), c(
    apply(d_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[2]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[2]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )

  plot(d_wv_std[[3]],
    type = "l", xlab = "Nível da variância de ondaletas", ylab = "Desvio padrão \n (n=2048)", ylim = c(0, 0.030),
    cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.sub = 1.5,
    lwd = 1.5
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col1, border = "NA"
  )
  polygon(c(1:8, 8:1), c(
    apply(d_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    rev(apply(d_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }))
  ),
  col = col2, border = "NA"
  )
  lines(d_wv_std[[3]],
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sb_8_std[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "dotted",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[2]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  lines(
    apply(d_wv_sd_multitaper[[3]], 2, function(x) {
      quantile(x)[[4]]
    }),
    lty = "longdash",
    lwd = 1.5
  )
  legend("topleft",
    legend = c("Multitaper", "Bootstrap"),
    fill = c(col1, col2), cex = 1.5
  )
  dev.off()
}

save.image(file.path(workspace_dir, "4_std_error.RData"))
