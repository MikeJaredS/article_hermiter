# Load libraries
options(dplyr.summarise.inform = F)
library(hermiter)
library(tdigest)
library(microbenchmark)
library(dplyr)
library(randtoolbox)
library(ggplot2)

# Set random seed for reproducibility
set.seed(20)

# Begin: R implementation of count matrix algorithm of Xiao, Wei. "Novel online
# algorithms for nonparametric correlations with application to analyze sensor
# data." 2019 IEEE International Conference on Big Data (Big Data). IEEE,
# as an S3 class. The implementation below follows 
# https://github.com/wxiao0421/onlineNPCORR/ reasonably closely in parts.

count_matrix_calculator <-
  function(cut_points_inpt, normalize = FALSE) {
    this <- list(
      cut_points = cut_points_inpt,
      count_matrix = matrix(
        rep(0, cut_points_inpt ^ 2),
        nrow = cut_points_inpt,
        ncol = cut_points_inpt,
        byrow = TRUE
      ),
      n_row = rep(0, cut_points_inpt),
      n_col = rep(0, cut_points_inpt),
      x_breaks = qnorm(p = seq(0, 1, length.out =
                                 (
                                   cut_points_inpt + 1
                                 )))[2:(cut_points_inpt)],
      y_breaks = qnorm(p = seq(0, 1, length.out =
                                 (
                                   cut_points_inpt + 1
                                 )))[2:(cut_points_inpt)],
      num_obs = 0,
      normalize_obs = normalize,
      running_mean_x = 0,
      running_mean_y = 0,
      running_variation_x = 0,
      running_variation_y = 0
    )
    class(this) <- append(class(this), "count_matrix_calculator")
    return(this)
  }

get_idx <- function(x, breaks) {
  if (breaks[length(breaks)] < x) {
    return(length(breaks) + 1)
  }
  return(which(breaks >= x)[1])
}

update_matrix <- function(this, x, y) {
  UseMethod("update_matrix", this)
}

update_matrix.count_matrix_calculator <- function(this, x) {
  this$num_obs <- this$num_obs + 1
  if (this$normalize_obs == TRUE) {
    prev_mean <- c(this$running_mean_x, this$running_mean_y)
    upd_mean <- (prev_mean * (this$num_obs - 1) + x) / this$num_obs
    this$running_mean_x  <- upd_mean[1]
    this$running_mean_y  <- upd_mean[2]
    if (this$num_obs < 2) {
      return(this)
    }
    upd_var <- c(this$running_variation_x,
                 this$running_variation_y) + (x - prev_mean) *
      (x - upd_mean)
    this$running_variation_x <- upd_var[1]
    this$running_variation_y <- upd_var[2]
    x <- (x - upd_mean) / sqrt(upd_var / (this$num_obs))
  }
  idx_row <- get_idx(x[1], this$x_breaks)
  idx_col <- get_idx(x[2], this$y_breaks)
  this$count_matrix[idx_row, idx_col] <-
    this$count_matrix[idx_row, idx_col] + 1
  this$n_row[idx_row] <- this$n_row[idx_row] + 1
  this$n_col[idx_col] <- this$n_col[idx_col] + 1
  return(this)
}

get_spearmans <- function(this) {
  UseMethod("get_spearmans", this)
}

get_spearmans.count_matrix_calculator <- function(this) {
  len_x_breaks <- length(this$x_breaks)
  len_y_breaks <- length(this$y_breaks)
  r_row <- rep(0, len_x_breaks + 1)
  r <- 0
  for (k in c(1:(len_x_breaks + 1))) {
    if (this$n_row[k] == 0) {
      r_row[k] <- r
    } else {
      r_row[k] <- ((r + 1) + (r + this$n_row[k])) / 2
      r <- r + this$n_row[k]
    }
  }
  r_col <- rep(0, len_y_breaks + 1)
  r <- 0
  for (k in c(1:(len_y_breaks + 1))) {
    if (this$n_col[k] == 0) {
      r_col[k] <- r
    } else {
      r_col[k] <- ((r + 1) + (r + this$n_col[k])) / 2
      r <- r + this$n_col[k]
    }
  }
  r_row <- r_row - (this$num_obs + 1) / 2
  r_col <- r_col - (this$num_obs + 1) / 2
  r_row <- r_row / sqrt(sum(this$n_row * r_row ^ 2))
  r_col <- r_col / sqrt(sum(this$n_col * r_col ^ 2))
  corr <- t(r_row) %*% this$count_matrix %*% r_col
  return(corr)
}

get_kendall <- function(this) {
  UseMethod("get_kendall", this)
}

get_kendall.count_matrix_calculator <- function(this) {
  len_n_row <- length(this$n_row)
  len_n_col <- length(this$n_col)
  count_mat_sum <-
    matrix(rep(0, len_n_row * len_n_col), len_n_row, len_n_col, byrow = TRUE)
  for (i in 2:len_n_row) {
    count_mat_sum[i, 2:len_n_col] <-
      cumsum(this$count_matrix[(i - 1), 1:(len_n_col - 1)])
  }
  for (i in 2:len_n_row) {
    count_mat_sum[i, ] <- count_mat_sum[i, ] + count_mat_sum[(i - 1), ]
  }
  concord_pairs <- sum(this$count_matrix * count_mat_sum)
  ties_in_x <- 0
  for (i in 1:len_n_row) {
    ties_in_x <- ties_in_x + (this$n_row[i] ^ 2 - 
                                sum(this$count_matrix[i, ] ^ 2)) / 2
  }
  ties_in_y <- 0
  for (j in 1:len_n_col) {
    ties_in_y <- ties_in_y + (this$n_col[j] ^ 2 - 
                                sum(this$count_matrix[, j] ^ 2)) / 2
  }
  ties_in_x_and_y <- sum(this$count_matrix * (this$count_matrix - 1)) / 2
  discord_pairs <- this$num_obs * (this$num_obs - 1) / 2 - concord_pairs -
    ties_in_x - ties_in_y - ties_in_x_and_y
  corr <- (concord_pairs - discord_pairs) / 
    sqrt((concord_pairs + discord_pairs + ties_in_x) * 
           (concord_pairs + discord_pairs + ties_in_y))
  return(corr)
}
# End: R implementation of count matrix algorithm of cite as an S3 class

# Begin: Univariate simulation study comparing hermiter and tdigest for quantile
# estimation
calculate_miae <- function(full_miae = FALSE) {
  distros_index <- c(1:5, 7:8, 11, 13:17, 21:28)
  numruns <- 1000
  num_obs_vec <- c(1e4, 1e5, 1e6, 1e7)
  if (full_miae == TRUE) {
    p <- randtoolbox::sobol(1000)
    norm_factor <- 1
  } else {
    p <- randtoolbox::sobol(1000) * 0.98 + 0.01
    norm_factor <- 0.98
  }
  distr_name_all <- c()
  num_obs_all <- c()
  mae_hermite_quant <- c()
  mae_t_digest_quant <- c()
  count <- 0
  for (num_obs in num_obs_vec) {
    for (current_distro_idx in seq_along(distros_index)) {
      dnum <- distros_index[current_distro_idx]
      print(current_distro_idx)
      distr_name <- benchden::berdev(dnum)$name
      r_func <-
        function(core_obs) {
          benchden::rberdev(n = core_obs, dnum = dnum)
        }
      q_func <-
        function(p_est) {
          benchden::qberdev(p_est, dnum = dnum)
        }
      p_func <- function(x) {
        benchden::pberdev(x, dnum = dnum)
      }
      res_hermite_quant <- rep(0, numruns)
      res_t_digest_quant <- rep(0, numruns)
      for (run in c(1:numruns)) {
        obs <- r_func(num_obs)
        h_est <- hermite_estimator(N = 30, standardize = T) %>%
          update_batch(obs)
        td <- tdigest(obs)
        q_est_hermite <- h_est %>% quant(p)
        q_est_t_digest <- quantile(td, probs = p)
        true_quant <- q_func(p)
        res_hermite_quant[run] <- norm_factor *
          mean(abs(q_est_hermite - true_quant))
        res_t_digest_quant[run] <- norm_factor *
          mean(abs(q_est_t_digest - true_quant))
      }
      mae_herm_quant <- mean(res_hermite_quant)
      mae_t_dig_quant <- mean(res_t_digest_quant)
      count <- count + 1
      distr_name_all[count] <- distr_name
      num_obs_all[count] <- num_obs
      mae_hermite_quant[count] <- mae_herm_quant
      mae_t_digest_quant[count] <- mae_t_dig_quant
    }
  }
  result <- data.frame(
    distribution_name = distr_name_all,
    num_obs =
      num_obs_all,
    mae_hermite_quant,
    mae_t_digest_quant
  )
  result <- result %>% mutate(hermite_better_quant =
                                ifelse(mae_hermite_quant < mae_t_digest_quant, 
                                       1, 0))
  univariate_quantile_results <- result %>%
    group_by(num_obs) %>%
    summarise(
      num_herm_better = sum(hermite_better_quant),
      total_distros = n(),
      perc_herm_better =
        sum(hermite_better_quant) / n(),
      p_val =
        wilcox.test(
          mae_hermite_quant,
          mae_t_digest_quant,
          paired = TRUE,
          alternative = "less",
          mu = 0
        )$p.value
    )
  return(univariate_quantile_results)
}
univar_quant_results_partial <-
  calculate_miae(full_miae = FALSE)
univar_quant_results_full <- 
  calculate_miae(full_miae = TRUE)
print(univar_quant_results_partial)
print(univar_quant_results_full)

# End: Univariate simulation study comparing hermiter and tdigest for quantile
# estimation

# Begin: Benchmark hermiter vs tdigest, updating with 1e6 observations
h_est_10 <- hermite_estimator(N = 10, standardize = T)
h_est_20 <- hermite_estimator(N = 20, standardize = T)
h_est_30 <- hermite_estimator(N = 30, standardize = T)
obs <- rnorm(1e6)
bench_res <- microbenchmark::microbenchmark(
  t_digest = td <- tdigest(obs),
  hermite_N_10 = update_batch(h_est_10, obs),
  hermite_N_20 = update_batch(h_est_20, obs),
  hermite_N_30 = update_batch(h_est_30, obs),
  times = 100
)
autoplot(bench_res, log = TRUE)
bench_res
# End: Benchmark hermiter vs tdigest, updating with 1e6 observations

# Begin: Benchmark hermiter vs tdigest, quantile estimation
obs <- rnorm(1e6)
td <- tdigest(obs)
h_est <-
  hermite_estimator(N = 30, standardize = T)  %>% update_batch(obs)
p_1 <- 0.5
p_100 <- seq(0.01, 1, 0.01)
p_10000 <- seq(0.0001, 1, 0.0001)
p_100000 <- seq(0.00001, 1, 0.00001)
bench_res <- microbenchmark::microbenchmark(
  hermite_1_quantile =  quant(h_est, p = p_1),
  tdigest_1_quantile = quantile(td, probs = p_1),
  hermite_100_quantiles = quant(h_est, p = p_100),
  tdigest_100_quantiles = quantile(td, probs = p_100),
  hermite_10_000_quantiles = quant(h_est, p = p_10000),
  tdigest_10_000_quantiles = quantile(td, probs = p_10000),
  hermite_100_000_quantiles = quant(h_est, p = p_100000),
  tdigest_100_000_quantiles = quantile(td, probs = p_100000),
  times = 1e2
)
autoplot(bench_res,
         log = T,
         xlab = "Algorithm",
         ylab = "Time (millis)")
bench_res
# End: Benchmark hermiter vs tdigest, quantile estimation

# Begin: Bivariate simulation study comparing hermiter and count matrix
# algorithms for estimation of Spearman's Rho and Kendall Tau coefficients.
rho_inpt <- c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)
num_obs_inpt <- c(1e4, 5e4, 1e5)
runs <- 1000
mae_matrix_kendall <- c()
mae_hermite_kendall <- c()
mae_matrix_spear <- c()
mae_hermite_spear <- c()
rho_vec <- c()
num_obs_vec <- c()
sig_x <- 1
sig_y <- 1
for (num_obs in num_obs_inpt) {
  print(num_obs)
  for (rho in rho_inpt) {
    print(rho)
    for (j in 1:runs) {
      obs_mat <- mvtnorm::rmvnorm(
        n = num_obs,
        mean = rep(0, 2),
        sigma = matrix(
          c(sig_x ^ 2, rho * sig_x * sig_y,
            rho * sig_x * sig_y, sig_y ^ 2),
          nrow = 2,
          ncol = 2,
          byrow = T
        )
      )
      matrix_est_c30 <-
        count_matrix_calculator(cut_points_inpt = 30,
                                normalize = F)
      matrix_est_c100 <-
        count_matrix_calculator(cut_points_inpt = 100,
                                normalize = F)
      hermite_est <-
        hermite_estimator_bivar(N = 30, standardize = F)
      for (i in seq_len(nrow(obs_mat))) {
        matrix_est_c100 <- matrix_est_c100 %>% update_matrix(obs_mat[i, ])
        matrix_est_c30 <-
          matrix_est_c30 %>% update_matrix(obs_mat[i, ])
        hermite_est <-
          hermite_est %>% update_sequential(obs_mat[i, ])
      }
      kendall_est_matrix <-
        matrix_est_c100 %>% get_kendall()
      kendall_est_hermite <- hermite_est %>% kendall()
      spear_est_matrix <-
        matrix_est_c30 %>% get_spearmans()
      spear_est_hermite <- hermite_est %>% spearmans()
      kendall_true <- 2 / pi * asin(rho)
      spear_true <- cor(obs_mat, method = "spearman")[1, 2]
      mae_matrix_kendall <- append(mae_matrix_kendall,
                                   abs(kendall_est_matrix - kendall_true))
      mae_hermite_kendall <- append(mae_hermite_kendall,
                                    abs(kendall_est_hermite - kendall_true))
      mae_matrix_spear <- append(mae_matrix_spear,
                                 abs(spear_est_matrix - spear_true))
      mae_hermite_spear <- append(mae_hermite_spear,
                                  abs(spear_est_hermite - spear_true))
      num_obs_vec <- append(num_obs_vec, num_obs)
      rho_vec <- append(rho_vec, rho)
    }
  }
}
result <-
  data.frame(
    num_obs_vec,
    rho_vec,
    mae_matrix_kendall,
    mae_hermite_kendall,
    mae_matrix_spear,
    mae_hermite_spear
  )
summary_by_rho_and_num_obs <-
  result %>% 
  group_by(num_obs_vec, rho_vec) %>%
  summarise(
    mae_matrix_kendall = mean(mae_matrix_kendall),
    mae_hermite_kendall =
      mean(mae_hermite_kendall),
    mae_matrix_spear = mean(mae_matrix_spear),
    mae_hermite_spear = mean(mae_hermite_spear)
  )
summary_num_obs_kendall <- summary_by_rho_and_num_obs %>%
  group_by(num_obs_vec) %>% 
  summarise(
    mae_matrix_kendall_avg =
      mean(mae_matrix_kendall) * 100,
    mae_hermite_kendall_avg =
      mean(mae_hermite_kendall) * 100,
    sd_matrix_kendall =
      sd(mae_matrix_kendall) * 100,
    sd_hermite_kendall = sd(mae_hermite_kendall) * 100
  )
summary_num_obs_spear  <- summary_by_rho_and_num_obs %>%
  group_by(num_obs_vec) %>%
  summarise(
    mae_matrix_spear_avg =
      mean(mae_matrix_spear) * 100,
    mae_hermite_spear_avg = mean(mae_hermite_spear) * 100,
    sd_matrix_spear = sd(mae_matrix_spear) * 100,
    sd_hermite_spear =
      sd(mae_hermite_spear) * 100
  )
print(summary_num_obs_kendall)
print(summary_num_obs_spear)

# End: bivariate simulation study comparing hermiter and count matrix
# algorithms for estimation of Spearman's Rho and Kendall Tau coefficients.
