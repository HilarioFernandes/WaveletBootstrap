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

# <- SET THIS before running
base_path <- "C:/Users/Hilar/Projects/WaveletBootstrap"
workspace_dir <- file.path(base_path, "src", "WorkspaceData")
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)

source(file.path(base_path, "src", "1_Simulation_functions.R"))
source(file.path(base_path, "src", "2_Bootstrap_methods.R"))

# --- Testing Mode ---
test_mode <- TRUE
# --------------------

################################################################################

# Study for distribution/variance estimation

set.seed(42)

b <- if (test_mode) 5 else 100
iterations <- if (test_mode) 2 else 100

a_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
b_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
c_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)
d_wv <- list("128" = NULL, "512" = NULL, "2048" = NULL)

a_wv_wb_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_bw_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_wb_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

a_wv_bw_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_wb_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_bw_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_wb_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

b_wv_bw_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_wb_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_bw_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_wb_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

c_wv_bw_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_wb_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_bw_nbb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_wb_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)

d_wv_bw_sb <- list(
  "128" = vector("list", length = iterations),
  "512" = vector("list", length = iterations),
  "2048" = vector("list", length = iterations)
)


# Design of the simulation study:
# 100 iterations × 100 bootstrap reps × 3 sample sizes × 4 models.

for (iter in 1:iterations) {
  print(iter)

  for (i in 1:3) {

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

    f_nbb <- function(n) floor(4 * log2(n))
    f_sb <- function(n) 1 / (4 * log2(n))

    a_wv_wb_nbb[[i]][[iter]] <-
      bootstrap_wavelet(ya, "wb", TRUE, "NBB", f_nbb, b)
    a_wv_bw_nbb[[i]][[iter]] <-
      bootstrap_wavelet(ya, "bw", TRUE, "NBB", f_nbb, b)
    a_wv_wb_sb[[i]][[iter]] <-
      bootstrap_wavelet(ya, "wb", TRUE, "SB", f_sb, b)
    a_wv_bw_sb[[i]][[iter]] <-
      bootstrap_wavelet(ya, "bw", TRUE, "SB", f_sb, b)

    b_wv_wb_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yb, "wb", TRUE, "NBB", f_nbb, b)
    b_wv_bw_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yb, "bw", TRUE, "NBB", f_nbb, b)
    b_wv_wb_sb[[i]][[iter]] <-
      bootstrap_wavelet(yb, "wb", TRUE, "SB", f_sb, b)
    b_wv_bw_sb[[i]][[iter]] <-
      bootstrap_wavelet(yb, "bw", TRUE, "SB", f_sb, b)

    c_wv_wb_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yc, "wb", TRUE, "NBB", f_nbb, b)
    c_wv_bw_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yc, "bw", TRUE, "NBB", f_nbb, b)
    c_wv_wb_sb[[i]][[iter]] <-
      bootstrap_wavelet(yc, "wb", TRUE, "SB", f_sb, b)
    c_wv_bw_sb[[i]][[iter]] <-
      bootstrap_wavelet(yc, "bw", TRUE, "SB", f_sb, b)

    d_wv_wb_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yd, "wb", TRUE, "NBB", f_nbb, b)
    d_wv_bw_nbb[[i]][[iter]] <-
      bootstrap_wavelet(yd, "bw", TRUE, "NBB", f_nbb, b)
    d_wv_wb_sb[[i]][[iter]] <-
      bootstrap_wavelet(yd, "wb", TRUE, "SB", f_sb, b)
    d_wv_bw_sb[[i]][[iter]] <-
      bootstrap_wavelet(yd, "bw", TRUE, "SB", f_sb, b)
  }
}

################################################################################

# Distribution comparisons

#' Calculate maximum absolute difference between empirical cdfs
#'
#' @param list1 First nested list of wavelet variance estimates
#' @param list2 Second nested list of bootstrap wavelet variance estimates
#' @return A list of matrices (one per sample size) containing the max abs diff
#'   for each iteration and level
max_abs_diff_cdfs <- function(list1, list2) {
  # We create a matrix in which rows correspond to iterations and columns
  # correspond to scales. Each element of this matrix is the maximum absolute
  # difference between the bootstrap empirical cdfs

  list_temp <- list(
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (128 - 1) / (8 - 1)))),
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (512 - 1) / (8 - 1)))),
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (2048 - 1) / (8 - 1))))
  )


  for (i in 1:3) {
    for (iter in 1:iterations) {
      # we check the maximum level for which there is at least one bootstrap
      # estimate on that particular iteration
      n_levels <- max(which(apply(list2[[i]][[iter]], 2, function(x) {
        sum(1 - is.na(x))
      }) >= 2))

      for (j in 1:n_levels) {
        # We will then fill the matrices
        list_temp[[i]][iter, j] <- distance_cdfs( # nolint: object_usage_linter
          list2[[i]][[iter]][, j], list1[[i]][, j], 100
        )
      }
    }
  }

  list_temp
}

# comparing true values with (bootstrap then wavelet transformation)

# a_wv vs a_wv_bw_nbb
a_nbb_bw_dist <- max_abs_diff_cdfs(a_wv, a_wv_bw_nbb)

# a_wv vs a_wv_bw_sb
a_sb_bw_dist <- max_abs_diff_cdfs(a_wv, a_wv_bw_sb)

# b_wv vs b_wv_bw_nbb
b_nbb_bw_dist <- max_abs_diff_cdfs(b_wv, b_wv_bw_nbb)

# b_wv vs b_wv_bw_sb
b_sb_bw_dist <- max_abs_diff_cdfs(b_wv, b_wv_bw_sb)

# a_wv vs c_wv_bw_nbb
c_nbb_bw_dist <- max_abs_diff_cdfs(c_wv, c_wv_bw_nbb)

# c_wv vs c_wv_bw_sb
c_sb_bw_dist <- max_abs_diff_cdfs(c_wv, c_wv_bw_sb)

# d_wv vs d_wv_bw_nbb
d_nbb_bw_dist <- max_abs_diff_cdfs(d_wv, d_wv_bw_nbb)

# d_wv vs d_wv_bw_sb
d_sb_bw_dist <- max_abs_diff_cdfs(d_wv, d_wv_bw_sb)


# comparing true values with (wavelet transformation then bootstrap)

# a_wv vs a_wv_wb_nbb
a_nbb_wb_dist <- max_abs_diff_cdfs(a_wv, a_wv_wb_nbb)

# a_wv vs a_wv_wb_sb
a_sb_wb_dist <- max_abs_diff_cdfs(a_wv, a_wv_wb_sb)

# b_wv vs b_wv_wb_nbb
b_nbb_wb_dist <- max_abs_diff_cdfs(b_wv, b_wv_wb_nbb)

# b_wv vs b_wv_wb_sb
b_sb_wb_dist <- max_abs_diff_cdfs(b_wv, b_wv_wb_sb)

# a_wv vs c_wv_wb_nbb
c_nbb_wb_dist <- max_abs_diff_cdfs(c_wv, c_wv_wb_nbb)

# c_wv vs c_wv_wb_sb
c_sb_wb_dist <- max_abs_diff_cdfs(c_wv, c_wv_wb_sb)

# d_wv vs d_wv_wb_nbb
d_nbb_wb_dist <- max_abs_diff_cdfs(d_wv, d_wv_wb_nbb)

# d_wv vs d_wv_wb_sb
d_sb_wb_dist <- max_abs_diff_cdfs(d_wv, d_wv_wb_sb)


################################################################################

# Variance comparisons

#' Calculate ratio of variances
#'
#' @param list1 First nested list of wavelet variance estimates
#' @param list2 Second nested list of bootstrap wavelet variance estimates
#' @return A list of matrices with minimum of var ratio and its inverse
ratio_vars <- function(list1, list2) {
  # We create a matrix in which rows correspond to iterations and columns
  # correspond to scales. Each element of this matrix is the minimum between
  # the ratio of the corresponding bootstrap variance estimates and the
  # inverse ratio

  list_temp <- list(
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (128 - 1) / (8 - 1)))),
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (512 - 1) / (8 - 1)))),
    matrix(NA, nrow = iterations, ncol = floor(log2(1 + (2048 - 1) / (8 - 1))))
  )


  for (i in 1:3) {
    for (iter in 1:iterations) {
      # we check the maximum level for which there is at least one bootstrap
      # estimate on that particular iteration
      n_levels <- max(which(apply(list2[[i]][[iter]], 2, function(x) {
        sum(1 - is.na(x))
      }) >= 2))

      for (j in 1:n_levels) {
        ratio <- var(list2[[i]][[iter]][, j]) / var(list1[[i]][, j])

        # We will then fill the matrices
        list_temp[[i]][iter, j] <- min(ratio, 1 / ratio)
      }
    }
  }

  list_temp
}

# comparing true values with (bootstrap then wavelet transformation)

# a_wv vs a_wv_bw_nbb
a_nbb_bw_var <- ratio_vars(a_wv, a_wv_bw_nbb)

# a_wv vs a_wv_bw_sb
a_sb_bw_var <- ratio_vars(a_wv, a_wv_bw_sb)

# b_wv vs b_wv_bw_nbb
b_nbb_bw_var <- ratio_vars(b_wv, b_wv_bw_nbb)

# b_wv vs b_wv_bw_sb
b_sb_bw_var <- ratio_vars(b_wv, b_wv_bw_sb)

# a_wv vs c_wv_bw_nbb
c_nbb_bw_var <- ratio_vars(c_wv, c_wv_bw_nbb)

# c_wv vs c_wv_bw_sb
c_sb_bw_var <- ratio_vars(c_wv, c_wv_bw_sb)

# d_wv vs d_wv_bw_nbb
d_nbb_bw_var <- ratio_vars(d_wv, d_wv_bw_nbb)

# d_wv vs d_wv_bw_sb
d_sb_bw_var <- ratio_vars(d_wv, d_wv_bw_sb)


# comparing true values with (wavelet transformation then bootstrap)

# a_wv vs a_wv_wb_nbb
a_nbb_wb_var <- ratio_vars(a_wv, a_wv_wb_nbb)

# a_wv vs a_wv_wb_sb
a_sb_wb_var <- ratio_vars(a_wv, a_wv_wb_sb)

# b_wv vs b_wv_wb_nbb
b_nbb_wb_var <- ratio_vars(b_wv, b_wv_wb_nbb)

# b_wv vs b_wv_wb_sb
b_sb_wb_var <- ratio_vars(b_wv, b_wv_wb_sb)

# a_wv vs c_wv_wb_nbb
c_nbb_wb_var <- ratio_vars(c_wv, c_wv_wb_nbb)

# c_wv vs c_wv_wb_sb
c_sb_wb_var <- ratio_vars(c_wv, c_wv_wb_sb)

# d_wv vs d_wv_wb_nbb
d_nbb_wb_var <- ratio_vars(d_wv, d_wv_wb_nbb)

# d_wv vs d_wv_wb_sb
d_sb_wb_var <- ratio_vars(d_wv, d_wv_wb_sb)

################################################################################

# Summarizing the results (distribution estimation)

#' Summarize distribution study results
#'
#' @param list1 A list of matrices generated by max_abs_diff_cdfs
#' @return A matrix of summarized quantiles
summary_study_dist <- function(list1) {
  temp <- matrix(NA, nrow = 15, ncol = 8)

  for (i in 1:3) {
    n_levels <- max(which(apply(list1[[i]], 2, function(x) {
      sum(1 - is.na(x))
    }) >= 2))

    for (j in 1:n_levels) {
      temp[(5 * (i - 1) + 1):(5 * i), j] <- quantile(
        list1[[i]][, j], na.rm = TRUE
      )
    }
  }

  round(temp, digits = 2)
}

# bootstrap then wavelet transform

a_nbb_bw_dist_summary <- summary_study_dist(a_nbb_bw_dist)

a_sb_bw_dist_summary <- summary_study_dist(a_sb_bw_dist)

b_nbb_bw_dist_summary <- summary_study_dist(b_nbb_bw_dist)

b_sb_bw_dist_summary <- summary_study_dist(b_sb_bw_dist)

c_nbb_bw_dist_summary <- summary_study_dist(c_nbb_bw_dist)

c_sb_bw_dist_summary <- summary_study_dist(c_sb_bw_dist)

d_nbb_bw_dist_summary <- summary_study_dist(d_nbb_bw_dist)

d_sb_bw_dist_summary <- summary_study_dist(d_sb_bw_dist)

# wavelet transform then bootstrap

a_nbb_wb_dist_summary <- summary_study_dist(a_nbb_wb_dist)

a_sb_wb_dist_summary <- summary_study_dist(a_sb_wb_dist)

b_nbb_wb_dist_summary <- summary_study_dist(b_nbb_wb_dist)

b_sb_wb_dist_summary <- summary_study_dist(b_sb_wb_dist)

c_nbb_wb_dist_summary <- summary_study_dist(c_nbb_wb_dist)

c_sb_wb_dist_summary <- summary_study_dist(c_sb_wb_dist)

d_nbb_wb_dist_summary <- summary_study_dist(d_nbb_wb_dist)

d_sb_wb_dist_summary <- summary_study_dist(d_sb_wb_dist)


################################################################################


# Summarizing the results (variance estimation)

#' Summarize variance study results
#'
#' @param list1 A list of matrices generated by ratio_vars
#' @return A matrix of summarized quantiles
summary_study_var <- function(list1) {
  temp <- matrix(NA, nrow = 15, ncol = 8)

  for (i in 1:3) {
    n_levels <- max(which(apply(list1[[i]], 2, function(x) {
      sum(1 - is.na(x))
    }) >= 2))

    for (j in 1:n_levels) {
      temp[(5 * (i - 1) + 1):(5 * i), j] <- quantile(
        list1[[i]][, j], na.rm = TRUE
      )
    }
  }

  round(temp, digits = 2)
}


# bootstrap then wavelet transform

a_nbb_bw_var_summary <- summary_study_var(a_nbb_bw_var)

a_sb_bw_var_summary <- summary_study_var(a_sb_bw_var)

b_nbb_bw_var_summary <- summary_study_var(b_nbb_bw_var)

b_sb_bw_var_summary <- summary_study_var(b_sb_bw_var)

c_nbb_bw_var_summary <- summary_study_var(c_nbb_bw_var)

c_sb_bw_var_summary <- summary_study_var(c_sb_bw_var)

d_nbb_bw_var_summary <- summary_study_var(d_nbb_bw_var)

d_sb_bw_var_summary <- summary_study_var(d_sb_bw_var)

# wavelet transform then bootstrap

a_nbb_wb_var_summary <- summary_study_var(a_nbb_wb_var)

a_sb_wb_var_summary <- summary_study_var(a_sb_wb_var)

b_nbb_wb_var_summary <- summary_study_var(b_nbb_wb_var)

b_sb_wb_var_summary <- summary_study_var(b_sb_wb_var)

c_nbb_wb_var_summary <- summary_study_var(c_nbb_wb_var)

c_sb_wb_var_summary <- summary_study_var(c_sb_wb_var)

d_nbb_wb_var_summary <- summary_study_var(d_nbb_wb_var)

d_sb_wb_var_summary <- summary_study_var(d_sb_wb_var)

################################################################################

#' Generate LaTeX output for summary matrices
#'
#' @param matrix1 First summary matrix
#' @param matrix2 Second summary matrix
#' @return Prints formatted LaTeX table rows
latex_output <- function(matrix1, matrix2) {
  for (i in 1:15) {
    temp1 <- ""
    temp2 <- ""

    if ((i - 1) %% 5 == 2) {
      temp1 <- 2^((2 * (1 + floor((i - 1) / 5)) - 1) + 6)
    }

    temp1 <- paste0(temp1, " & ", ((i - 1) %% 5) * 25, "\\% ")
    temp2 <- paste0(temp2, " & ")

    for (j in 1:(2 * (2 + floor((i - 1) / 5)))) {
      temp1 <- paste0(temp1, " & ", matrix1[i, j])

      temp2 <- paste0(temp2, "& (", matrix2[i, j], ") ")
    }

    temp1 <- paste0(temp1, " \\\\ ")
    temp2 <- paste0(temp2, " \\\\ ")


    if ((i - 1) %% 5 == 4) {
      temp2 <- paste0(temp2, " \\hline")
    }

    cat(paste0(temp1, " \n"))
    cat(paste0(temp2, " \n"))
  }
}

latex_output(a_nbb_bw_dist_summary, a_nbb_wb_dist_summary)
latex_output(a_sb_bw_dist_summary, a_sb_wb_dist_summary)

latex_output(b_nbb_bw_dist_summary, b_nbb_wb_dist_summary)
latex_output(b_sb_bw_dist_summary, b_sb_wb_dist_summary)

latex_output(c_nbb_bw_dist_summary, c_nbb_wb_dist_summary)
latex_output(c_sb_bw_dist_summary, c_sb_wb_dist_summary)

latex_output(d_nbb_bw_dist_summary, d_nbb_wb_dist_summary)
latex_output(d_sb_bw_dist_summary, d_sb_wb_dist_summary)


latex_output(a_nbb_bw_var_summary, a_nbb_wb_var_summary)
latex_output(a_sb_bw_var_summary, a_sb_wb_var_summary)

latex_output(b_nbb_bw_var_summary, b_nbb_wb_var_summary)
latex_output(b_sb_bw_var_summary, b_sb_wb_var_summary)

latex_output(c_nbb_bw_var_summary, c_nbb_wb_var_summary)
latex_output(c_sb_bw_var_summary, c_sb_wb_var_summary)

latex_output(d_nbb_bw_var_summary, d_nbb_wb_var_summary)
latex_output(d_sb_bw_var_summary, d_sb_wb_var_summary)

save.image(file.path(workspace_dir, "3_Bootstrap_Wavelet_order.RData"))
