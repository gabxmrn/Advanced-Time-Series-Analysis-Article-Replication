bs_main <- function(a, b, y, x, u, prior_specifications) {

  D <- u %*% t(u) / ncol(y)
  inv_A <- solve(a)

  kappa <- 0.5

  # Compute omega
  omega <- inv_A %*% D %*% t(inv_A)

  # Compute somega
  ar_omega <- aromega_computation(y, p = 2)





  # test <- log_likelihood(a, y, kappa, omega, somega)

  return(ar_omega)
}

mh <- function() {

}

posterior_a <- function(a, prior_spec, weight_matrix) {

  prior <- sum_prior_a(a, prior_spec, weight_matrix)
  likelihood <- 1

  posterior <- prior + likelihood

  return(posterior)
}

log_likelihood <- function(A, Y, kappa, omega, ar_omega) {

  t <- ncol(Y)

  num <- t * 0.5 * log(det(t(A) %*% omega %*% A))
  tau <- diag(kappa) * diag(t(A) %*% ar_omega %*% A)

  return(tau)
}

sum_prior_a <- function(A, prior_spec, weight_matrix) {

  # Get model parameters - a_temp_gdp and a_prec_gdp
  mode_tg_pg <- 0
  scale_tg_pg <- 0.01
  dof_tg_pg <- 3

  # Loop on countries
  countries <- as.character(row.names(prior_spec))

  # Create temporary vector for result
  temp_prior <- c()

  for (country in countries) {

    # Get elements of matrix a
    name_temp <- paste0(country, ".temp")
    name_prec <- paste0(country, ".prec")
    name_gdp <- paste0(country, ".gdp")

    a_tg <- A[name_temp, name_gdp]
    a_pg <- A[name_prec, name_gdp]
    a_gt <- A[name_gdp, name_temp]
    a_gp <- A[name_gdp, name_prec]

    # Compute h_d
    det_a <- 1 + a_pg * a_gp - a_tg * a_gt
    h_d_gt <- a_gt / det_a
    h_d_gp <- a_gp / det_a

    # Compute h_f
    w <- rowSums(weight_matrix[country, ])

    h_f_gt <- h_d_gt * w
    h_f_gp <- h_d_gp * w

    # Get model parameters - a_gdp_temp and a_gdp_prec
    mode_gt_gp <- prior_spec[country, 1]
    scale_gt_gp <- prior_spec[country, 2]
    dof_gt_gp <- prior_spec[country, 3]

    # Get model parameters - h_f
    mode_temp <- prior_spec[, "mode"]
    row_temp <- as.numeric(weight_matrix[country, ])

    mode_h_f <- crossprod(row_temp, mode_temp)
    scale_h_f <- 0.01
    dof_h_f <- 3

    # Compute prior for country
    prior <- log(prior_student(a_tg, mode_tg_pg, dof_tg_pg, scale_tg_pg)) +
      log(prior_student(a_pg, mode_tg_pg, dof_tg_pg, scale_tg_pg)) +
      log(prior_student(a_gt, mode_gt_gp, dof_gt_gp, scale_gt_gp)) +
      log(prior_student(a_gp, mode_gt_gp, dof_gt_gp, scale_gt_gp)) +
      log(prior_student(h_d_gt, mode_gt_gp, dof_gt_gp, scale_gt_gp)) +
      log(prior_student(h_d_gp, mode_gt_gp, dof_gt_gp, scale_gt_gp)) +
      log(prior_student(h_f_gt, mode_h_f, dof_h_f, scale_h_f)) +
      log(prior_student(h_f_gp, mode_h_f, dof_h_f, scale_h_f))

    # Stock results in temporary df
    temp_prior <- c(temp_prior, prior)

  }

  # Create result dataframe
  result <- as.data.frame(temp_prior)
  rownames(result) <- countries

  print(min(temp_prior))
  print(max(temp_prior))

  return(result)
}

prior_student <- function(x, mu, nu, sigma) {

  library(LaplacesDemon)
  prior <- dst(x = x, mu = mu, sigma = sigma, nu = nu, log = FALSE)

  return(prior)
}

aromega_computation <- function(y, p) {

  # Result matrix
  error_df <- c()

  n <- ncol(y)

  for (i in seq_len(nrow(y))) {

    y_temp <- as.matrix(y[i, (1 + p):n])
    y_t <- t(y[i, ])

    x_temp <- do.call(cbind, lapply(1:p, function(lag) y_t[(p + 1 - lag):(n - lag), ]))
    x_temp <- t(x_temp)
    x_temp <- as.matrix(rbind(x_temp, 1))

    # Compute beta
    beta <- y_temp %*% t(x_temp) %*% solve(x_temp %*% t(x_temp))

    # Compute error
    error <- y_temp - beta %*% x_temp

    # Stock error
    error_df <- rbind(error_df, error)
  }

  # Compute variance-covariance matrix
  var_cov <- error_df %*% t(error_df) / (n - p)

  return(var_cov)
}