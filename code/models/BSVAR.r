bs_main <- function(svar_model, prior_specifications) {

  # Get output matrix for SVARX model
  A <- svar_model[[1]]
  B <- svar_model[[2]]
  Y <- svar_model[[3]]
  X <- svar_model[[4]]
  U <- svar_model[[5]]

  inv_A <- solve(A) # Inverse matrix of A
  D <- U %*% t(U) / ncol(Y) # Variance-covariance matrix

  kappa <- 0.5 # Value set by article

  # Compute omega
  omega <- inv_A %*% D %*% t(inv_A)

  # Compute somega
  ar_omega <- aromega_computation(Y, p = 1) # warning p=2 chez B&H


  # Launch MH algorithm -> gives all matrix A and zeta for all draws

  # Compute all diagonal matrices D (Inverse gamma distribution needing kappa_star (consant) and zeta_star)

  # Compute all matrices B (normal with formulas for mean and var-cov)

  # Launch IRF estimations (needs horizon, all matrices A and B, interval and if it is cumulative)
  # Two sets of IRF: real GDP growth on precipitation or temperature shocks
  # Steps:
    # loop on number of draws
    # Compute H and Phi (B*H) because we are in the case of eq (11)
    # Compute IRF (see method in code)
    # Case if the IRF function is cumulative
    # Store IRF
    # Once the loop for the country has ended -> compute lower, upper bound (depending on CI) and median
    # Should output: for each country, for each time period: lower and upper bound + median

  # Plot IRF + make a map for better visualisation

  test <- log_likelihood(A, Y, kappa, omega, ar_omega)
  #test2 <- compute_zeta_star(A, X, Y, ar_omega, p = 2)

  return(test)
}

mh <- function(iter, burnin) {

  # Output vectors for matrices A and zeta

  # Compute posterior for initial matrix

  # Loop on number of iterations

    # Generate a new A matrix after a shock -> needs a shock function to estimate H_i

    # Compute the beta and zeta associated to the new matrix

    # Compute posterior of new matrix

    # Compute acceptance ratio

    # Acceptance/rejection step depending on threshold

      # Acceptance: keep proposal

      # Rejection: we don't keep the proposal and start again
    
    # If iter > burnin then stock A and zeta in output matrix

}

update_a <- function(A, mu, weight_matrix) {

  H <- solve(A)

  # To do :)
  H_F <- 1

  H_updated <- H + mu * H_F
  A_updated <- solve(H_updated)

  return(A_updated)
}

posterior_a <- function(A, Y, kappa, omega, ar_omega, prior_spec, weight_matrix) {

  prior <- sum_prior_a(A, prior_spec, weight_matrix)
  likelihood <- log_likelihood(A, Y, kappa, omega, ar_omega)

  posterior <- prior + likelihood

  return(posterior)
}

log_likelihood <- function(A, Y, kappa, omega, ar_omega) {

  t <- ncol(Y) # Number of observations
  A <- as.matrix(A) # Convert A from a df to a matrix

  # Compute numerator of log-likelihood
  num <- t * 0.5 * log(det(t(A) %*% omega %*% A))
  tau <- kappa * diag(t(A) %*% ar_omega %*% A)

  for (i in length(tau)) {
    num <- num + kappa * log(tau[i])
  }

  # Compute denominator of log-likelihood
  kappa_star <- kappa + t / 2
  #tau_star <- tau + zeta_star / 2
  #denom <-

  # Compute log-likelihood
  # like <- num - denom

  return(tau)
}

compute_zeta_star <- function(A, X, Y, ar_omega, p) {

  # Hyperparameters value
  lambda0 <- 0.5
  lambda1 <- 1
  lambda3 <- 100
  lambda4 <- 100

  # Element v1
  v1 <- matrix(data = (1:p), nrow = p, ncol = 1) ^ (-2 * lambda1)

  # Element v2
  v2 <- matrix(data = diag(solve(diag(diag(ar_omega)))), ncol = 1)

  # Element v3
  v1_kro_v2 <- kronecker(v1, v2)
  v3_temp <- c(lambda3 ^ 2, lambda4 ^ 2)

  v3 <- matrix(nrow = 0, ncol = 1)
  i <- 1

  while (i <= nrow(v1_kro_v2)) {
    temp <- v1_kro_v2[i:min(i + 5, nrow(v1_kro_v2)), , drop = FALSE]
    v3 <- rbind(v3, temp, matrix(v3_temp, nrow = 2, ncol = 1, byrow = TRUE))
    i <- i + 6
  }

  v3 <- lambda0 ^ 2 * v3

  # Compute matrices M and M^(-1)
  M <- diag(x = c(v3), nrow = nrow(v3), ncol = nrow(v3))
  inv_M <- solve(M)

  # Compute P, the Cholesky factor of M^(-1)
  P <- chol(inv_M)

  # Regression with Y_tild and X_tild
  temp_output <- c()
  c_temp <- 1
  c_temp2 <- 1

  for (i in seq_len(33)) {

    Y_temp <- as.matrix(Y[c_temp:(c_temp + 2), ])
    X_temp <- as.matrix(X[c_temp2:(c_temp2 + 7), ])

    P_i <- as.matrix(P[(8 * i - 7):(8 * i), (8 * i - 7):(8 * i)])

    for (j in seq_len(3)) {

      # Get correct element of A
      a_i_j <- as.matrix(A[3 * (i - 1) + j, (3 * i - 2):(3 * i)])

      # Update Y_tild
      Y_tild <- a_i_j %*% Y_temp
      v_temp <- matrix(0, 1, 8)
      Y_tild <- t(cbind(Y_tild, v_temp))

      # Update X_tild
      X_tild <- t(X_temp)
      X_tild <- rbind(X_tild, P_i)
      X_tild_t <- t(X_tild)

      # Compute beta
      beta <- solve(X_tild_t %*% X_tild) %*% X_tild_t %*% Y_tild

      # Compute residuals
      res <- Y_tild - X_tild %*% beta

      # Squared residuals
      res_sq <- t(res) %*% res

      # Store squared residuals

    }

    c_temp <- c_temp + 3
    c_temp2 <- c_temp2 + 8
  }

  # Compute variance-covariance matrix: zeta_star

  # Output - zeta_star
  return(res_sq)
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

  # Outout matrix
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