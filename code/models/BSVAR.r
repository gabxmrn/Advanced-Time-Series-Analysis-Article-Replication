bs_main <- function(svar_model, prior_specifications, weight_matrix) {

  #! Points d'attention
    #! Omega -> voir page 15/16 si bonne forme
    #! prior(A) -> somme des priors de chaque pays ?
    #! B does not depend on D -> sigma is the ar_omega
    #! draw de B -> normal univariate mais est-ce qu'on doit faire une multi variée ?

  # Output from SVARX model
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
  omega_test <- matrix(0, nrow(omega), ncol(omega))
  for (i in 1:(nrow(omega) %/% 3)) {
    omega_test[(3 * i - 2):(3 * i), (3 * i - 2):(3 * i)] <- omega[(3 * i - 2):(3 * i), (3 * i - 2):(3 * i)]
  }

  # Compute ar_omega
  ar_omega <- aromega_computation(Y, p = 1) # warning p=2 chez B&H

  # Algorithm parameters
  iter <- 10
  burnin <- 5

  # When drawing matrix A, positions to update
  position_update <- matrix(data = 0, nrow = 3, ncol = 3)
  position_update[3, 1] <- 1
  position_update[3, 2] <- 1
  position_update[1, 3] <- 1
  position_update[2, 3] <- 1

  # Draw matrices A
  temp_draw <- metropolis_hastings_algorithm(iter, burnin, A, X, Y, kappa, omega_test, ar_omega, prior_specifications, weight_matrix, position_update)
  A_draw <- temp_draw[[1]]
  zeta_draw <- temp_draw[[2]]

  # Update matrices A with information from H
  A_revised <- list()
  for (i in seq_along(A_draw)) {
    A_revised[[i]] <- revise_A(A_draw[[i]], mu = 1, weight_matrix)
  }

  # Draw matrices D
  D_draw <- draw_D(A_revised, zeta_draw, kappa, 27, ar_omega)

  # Draw matrices B
  B_draw <- draw_B(A_revised, X, Y, ar_omega, p = 2)

  test <- IRF_preprocessing(A_revised, B_draw, weight_matrix)

  return(test)
}

IRF_preprocessing <- function(A,B,wm) { 

  countries <- as.character(row.names(wm))

  for (country in countries) {

    A_i <- list()
    B_i <- list()

    c <- 1
    c2 <- 1

    for (i in seq_along(A)) {
      A_temp <- A[[i]]
      B_temp <- B[[i]]

      A_i[[i]] <- A_temp[c:(c+2), c:(c+2)]
      B_i[[i]] <- B_temp[c:(c+2), c2:(c2+7)]
        dimnames(A_i[[i]]) <- list(
    rownames(A_temp)[c:(c+2)],
    colnames(A_temp)[c:(c+2)]
  )
    }

    #Ligne pour lancer l'IRF (voir C)
    test <- IRF(A_i, B_i, 5, FALSE)

    c <- c + 3
    c2 <- c2 + 8

    return(A_i[[1]])

  }
}


IRF <- function(A_draw, B_draw, horizon, is_cumulative = FALSE) { # Launch IRF estimations (needs horizon, all matrices A and B, interval and if it is cumulative)
  # Initialization
  irf_results <- list()
  #Two sets of IRF: real GDP growth on precipitation or temperature shocks
  # Steps:
    #For each country, boucle sur le nombre de tirages
    
    # loop on number of draws
     for (i in seq_along(A_draw)) {
    A <- A_draw[[i]]
    B <- B_draw[[i]]
    # Compute H and Pi (B*H) because we are in the case of eq (11)
    H <- solve(A)
    Pi <- H %*% B
    colnames(Pi) <- colnames(B)
    # Compute IRF (see method in code): 
    # Initialization for every horizon: 
    # irf_temp <- matrix(0, nrow = nrow(A), ncol = horizon + 1)
    # shock <- diag(D)
  
    # irf_temp[, 1] <- H %*% shock
    # # Propagation sur tous les horizons
    # for (h in 1:horizon) {
    #   Pi_h <- shock %*% Pi 
    #   # irf_temp[, h + 1] <- t(Pi_h)
    #   # Pi <- Pi %*% B  #
    # }
    # # Case if the IRF function is cumulative
    # if (is_cumulative) {
    #   irf_temp <- t(apply(irf_temp, 1, cumsum))
    # }
    # # Store IRF
    #  irf_results[[i]] <- irf_temp
  }
    # Once the loop for the country has ended -> compute lower, upper bound (depending on CI) and median
    # Should output: for each country, for each time period: lower and upper bound + median

  # Plot IRF + make a map for better visualisation
#library (ggplot2)
#horizon <- 10
#irf_var1 <- sapply(results, , )

return(Pi)

}

draw_B <- function(A_revised, X, Y, ar_omega, p) {

  # Load library (function used: bdiag)
  library("Matrix")

  B_output <- list() # Output matrix

  for (i in seq_along(A_revised)) {

    # Get values from draw
    A <- A_revised[[i]]

    # Tilded regression
    regr <- compute_tilded_regr(A, X, Y, ar_omega, p = 2)

    # Get m_star and M_star
    m_star <- as.matrix(regr[[2]])
    M_star <- as.matrix(regr[[3]])

    #Initialize B matrix for draw i
    B_draw <- matrix(nrow = 0, ncol = 0)

    # Loop countries
    for (k in seq_len(33)) {

      B_temp <- c()

      # Loop on series
      for (j in seq_len(3)) {

        # Get parameters for normal distribution
        m_star_j <- m_star[[(3 * (k - 1) + j)]]
        M_star_j <- M_star[[(3 * (k - 1) + j)]]

        # Draw values from a normal univariate distribution
        B_row_j_draw <- normal(m_star_j, M_star_j)

        # Add to country values
        B_temp <- rbind(B_temp, B_row_j_draw)
      }

      B_draw <- bdiag(B_draw, B_temp)
    }

    # Output matrix B
    B_output[[i]] <- as.matrix(B_draw)
  }

  return(B_output)
}

normal <- function(mean, cov_var) {

  # Output vector
  output <- c()

  for (i in seq_along(mean)) {

    draw <- rnorm(1, mean = mean[i], sd = sqrt(cov_var[i, i]))

    output <- c(output, draw)
  }

  return(output)
}

draw_D <- function(A_draw, zeta_star_draw, kappa, nb_obs, ar_omega) {

  output_D <- list() # Output list

  # Compute kappa_star
  kappa_star <- kappa + nb_obs / 2

  # Loop on the draw
  for (i in seq_along(A_draw)) {

    # Get values from draw
    A <- A_draw[[i]]
    zeta_star <- zeta_star_draw[[i]]

    # Compute tau
    tau <- kappa * diag(t(A) %*% ar_omega %*% A)

    # Compute tau_star
    tau_star <- tau + zeta_star / 2

    # Draw diagonal values of D from inverse gamma
    draw_D <- inverse_gamma(kappa_star, tau_star)

    # Compute D
    D <- diag(draw_D)

    # Output
    output_D[[i]] <- D
  }

  return(output_D)
}

inverse_gamma <- function(kappa_star, tau_star) {

  # Output vector
  draws <- c()

  # Draw from gamma distribution
  for (i in seq_along(tau_star)) {

    gamma_draw <- rgamma(1, kappa_star, tau_star[[i]])
    inv_gamma_draw <- 1 / gamma_draw

    draws <- c(draws, inv_gamma_draw)
  }
  return(draws)
}

revise_A <- function(A, mu, weight_matrix) {

  # Load library (function used: bdiag)
  library("Matrix")

  H <- solve(A) # Inverse of A
  countries <- as.character(row.names(weight_matrix)) # List of countries

  h_f_i <- matrix(data = 0, nrow = 3, ncol = 3) # Temporary output
  H_F <- matrix(nrow = 0, ncol = 0) # Final matrix H_F

  # Loop on countries to compute h_f_i
  for (country in countries) {

    # Get elements of matrix a
    name_temp <- paste0(country, ".temp")
    name_prec <- paste0(country, ".prec")
    name_gdp <- paste0(country, ".gdp")

    a_tg <- - A[name_temp, name_gdp]
    a_pg <- - A[name_prec, name_gdp]
    a_gt <- - A[name_gdp, name_temp]
    a_gp <- - A[name_gdp, name_prec]

    # Compute h_d
    det_a <- 1 + a_pg * a_gp - a_tg * a_gt
    h_d_gt <- a_gt / det_a
    h_d_gp <- a_gp / det_a

    # Compute h_f
    w <- rowSums(weight_matrix[country, ])

    h_f_gt <- h_d_gt * w
    h_f_gp <- h_d_gp * w

    # Compute matrix H_F_i
    h_f_i[3, 1] <- h_f_gt
    h_f_i[3, 2] <- h_f_gp

    H_F <- bdiag(H_F, h_f_i)
  }

  H_updated <- as.matrix(H) + as.matrix(H_F)
  A_updated <- solve(H_updated)

  return(A_updated)
}

metropolis_hastings_algorithm <- function(iter, burnin, A_ini, X, Y, kappa, omega, ar_omega, prior_spec, weight_matrix, pU) {

  # Initialization - GVARX is the starting point
  A_old <- A_ini
  regr <- compute_tilded_regr(A_old, X, Y, ar_omega, p = 2)
  zeta_old <- regr[[1]]
  A_old_post <- posterior_A(A_old, X, Y, kappa, zeta_old, omega, ar_omega, prior_spec, weight_matrix)

  # Output vectors for matrices A and zeta
  output_A <- list()
  output_zeta <- list()

  # Loop on number of iterations
  for (draw in seq_len(iter)) {

    # Generate a new A matrix
    A_new <- draw_A(A_old, pU, scale = 1, weight_matrix)

    # Compute the beta and zeta associated to the new matrix
    regr <- compute_tilded_regr(A_new, X, Y, ar_omega, p = 2)
    zeta_new <- regr[[1]]

    # Compute posterior of new matrix
    A_new_post <- posterior_A(A_new, X, Y, kappa, zeta_new, omega, ar_omega, prior_spec, weight_matrix)

    # Iteration of MH algorithm
    acceptance <- exp(A_new_post - A_old_post) # Acceptance probability
    threshold <- runif(1, min = 0.0, max = 1.0) # Threshold

    # Check acceptance
    if (threshold <= acceptance) {

      A_old <- A_new
      zeta_old <- zeta_new
      A_old_post <- A_new_post

    }

    # If iter > burnin then stock A and zeta in output matrix
    if (draw > burnin) {

      output_A[[draw - burnin]] <- A_old
      output_zeta[[draw - burnin]] <- zeta_old

    }
  }

  return(list(output_A, output_zeta))
}

draw_A <- function(A, pU, scale, wm) {

  # Load library (function used: bdiag)
  library("Matrix")

  A_updated <- matrix(nrow = 0, ncol = 0) # Output matrix

  countries <- as.character(rownames(wm))
  c <- 1

  for (country in countries) {

    # Get matrix A for the country
    a_i <- A[c:(c + 2), c:(c + 2)]
    c <- c + 3

    # Iterate on rows
    for (i in seq_len(3)) {

      # Iterate on columns
      for (j in seq_len(3)) {

        # Update value if necessary
        if (pU[i, j] == 1) {
          a_i[i, j] <- a_i[i, j] +
            scale * (rnorm(1) / sqrt(0.5 * (rnorm(1)^2 + rnorm(1)^2)))
        }
      }
    }

    # Add matrix a_i to output matrix
    A_updated <- bdiag(A_updated, as.matrix(a_i))
  }

  # Output formatting
  A_updated <- as.data.frame(as.matrix(A_updated))

  count2 <- 1
  for (country in countries) {
    rownames(A_updated)[count2] <- paste0(country, ".temp")
    rownames(A_updated)[count2 + 1] <- paste0(country, ".prec")
    rownames(A_updated)[count2 + 2] <- paste0(country, ".gdp")

    colnames(A_updated)[count2] <- paste0(country, ".temp")
    colnames(A_updated)[count2 + 1] <- paste0(country, ".prec")
    colnames(A_updated)[count2 + 2] <- paste0(country, ".gdp")

    count2 <- count2 + 3
  }

  return(as.matrix(A_updated))
}

posterior_A <- function(A, X, Y, kappa, zeta_star, omega, ar_omega, prior_spec, weight_matrix) {

  # Compute prior (using log)
  prior_a_country <- sum_prior_a(A, prior_spec, weight_matrix)
  prior_a <- sum(prior_a_country)

  # Compute log_likelihood
  likelihood <- log_likelihood(A, X, Y, kappa, zeta_star, omega, ar_omega)

  # Compute posterior for matrix A
  posterior <- prior_a + likelihood

  return(posterior)
}

log_likelihood <- function(A, X, Y, kappa, zeta_star, omega, ar_omega) {

  t <- ncol(Y) # Number of observations
  A <- as.matrix(A) # Convert A from a df to a matrix

  # Compute numerator of log-likelihood
  num <- t * 0.5 * log(det(t(A) %*% omega %*% A)) + kappa
  tau <- kappa * diag(t(A) %*% ar_omega %*% A)

  for (i in length(tau)) {
    num <- num + kappa * log(tau[i])
  }

  # Compute denominator of log-likelihood
  kappa_star <- kappa + t / 2
  tau_star <- tau + zeta_star / 2
  denom <- 0
  for (i in length(tau)) {
    denom <- denom + kappa_star * log((2 / t) * tau_star[i])
  }

  # Compute log-likelihood
  like <- num - denom

  return(like)
}

compute_tilded_regr <- function(A, X, Y, ar_omega, p) {

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
  c_temp <- 1
  c_temp2 <- 1
  n <- 1

  # Empty result vectors
  zeta <- c()
  m_star <- list()
  M_star <- list()

  for (i in seq_len(33)) {

    Y_temp <- as.matrix(Y[c_temp:(c_temp + 2), ])
    X_temp <- as.matrix(X[c_temp2:(c_temp2 + 7), ])

    P_i <- as.matrix(P[(8 * i - 7):(8 * i), (8 * i - 7):(8 * i)])

    for (j in seq_len(3)) {

      # Get correct element of A
      a_i_j <- matrix(A[3 * (i - 1) + j, (3 * i - 2):(3 * i)], nrow = 1, ncol = 3)

      # Update Y_tild
      Y_tild <- as.numeric(a_i_j) %*% Y_temp
      v_temp <- matrix(0, 1, 8)
      Y_tild <- t(cbind(Y_tild, v_temp))

      # Update X_tild
      X_tild <- t(X_temp)
      X_tild <- rbind(X_tild, P_i)
      X_tild_t <- t(X_tild)

      # Compute beta
      beta <- solve(X_tild_t %*% X_tild) %*% X_tild_t %*% Y_tild
      m_star[[n]] <- beta

      # Compute M_star
      M_temp <- solve(X_tild_t %*% X_tild)
      M_star[[n]] <- M_temp

      # Compute residuals
      res <- Y_tild - X_tild %*% beta

      # Squared residuals
      res_sq <- t(res) %*% res

      # Store squared residuals
      zeta <- c(zeta, res_sq)
      
      n <- n + 1
    }

    c_temp <- c_temp + 3
    c_temp2 <- c_temp2 + 8
  }

  # Output - zeta_star
  return(list(zeta, m_star, M_star))
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

    a_tg <- - A[name_temp, name_gdp]
    a_pg <- - A[name_prec, name_gdp]
    a_gt <- - A[name_gdp, name_temp]
    a_gp <- - A[name_gdp, name_prec]

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