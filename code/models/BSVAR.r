bs_main <- function(svar_model, prior_specifications, weight_matrix) {
  # ==========================================
  # Purpose: launches the bayesian algorithm
  # Parameters:
  #           svar_model: results the VARX, starting point for the bayesian algo
  #           prior_specifications: informations about the prior of some elements of A
  #           weight_matrix: trading weight matrix between the countries
  # Returns: nothing, plots an IRF
  # ==========================================

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
  ar_omega <- aromega_computation(Y, p = 1) # warning p = 2 chez B&H

  # Algorithm parameters
  iter <- 5 #4000
  burnin <- 1 #800

  # When drawing matrix A, positions to update
  position_update <- matrix(data = 0, nrow = 3, ncol = 3)
  position_update[3, 1] <- 1
  position_update[3, 2] <- 1
  position_update[1, 3] <- 1
  position_update[2, 3] <- 1

  # Draw matrices A
  temp_draw <- metropolis_hastings_algorithm(iter, burnin, A, X, Y, kappa,
                                             omega_test, ar_omega,
                                             prior_specifications,
                                             weight_matrix, position_update)
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

  # Launch IRF: impact of a shock on temperature or precipitation on gdp
  irf <- IRF_processing(A_revised, B_draw, weight_matrix)
  IRF_plot(irf, 5, "temp")
  IRF_plot(irf, 5, "prec")

  return(A_revised)
}

IRF_plot <- function(irf, horizon, series_shock) {
    # ==========================================
  # Purpose: plots the statistics of the IRFs for each country
  # Parameters:
  #           irf: a list of irfs
  #           horizon: time frame
  #           series_shock: which series has been shocked
  # Returns: nothing, prints a grid with the IRF of each country
  # ==========================================

  library(ggplot2)
  library(gridExtra)
  library(grid)

  irf_for_plot <- IRF_stats(irf, horizon, series_shock)

  create_plot <- function(df, name) {
    ggplot(df, aes(x = 1:nrow(df))) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.4) +
      geom_line(aes(y = median), color = "red", size = 1) +
      labs(x = "Years", y = "") +
      ggtitle(name) +
      theme_minimal() +
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            panel.grid = element_blank())
  }

  # Plot for each series
  plot_list <- lapply(names(irf_for_plot), function(name) {
    create_plot(irf_for_plot[[name]], name)
  })

  # Format global plot
  grid.arrange(grobs = plot_list, ncol = 6, nrow = 6)
  if (series_shock == "temp") {
    grid.text("Impact of a temperature shock on GDP", x = 0.5, y = 1,
              gp = gpar(fontsize = 20, fontface = "bold"))
  } else {
    grid.text("Impact of a precipitation shock on GDP", x = 0.5, y = 1,
              gp = gpar(fontsize = 20, fontface = "bold"))
  }
  return(irf_for_plot)
}

IRF_stats <- function(irf, horizon, series_shock) {
  # ==========================================
  # Purpose: computes the median, lower and upper bound of a set of IRF
  # Parameters:
  #           irf: a list of irfs
  #           horizon: time frame
  #           series_shock: which series has been shocked
  # Returns: the stats of the list of irf implemented
  # ==========================================

  if (series_shock == "temp") {
    ncol <- 3
  } else if (series_shock == "prec") {
    ncol <- 6
  }

  stats_df <- list()

  # Loop on countries
  countries <- names(irf)
  for (country in countries) {

    # Get list of irf
    irf_list <- irf[[country]]

    # Empty dataframes
    df_temp <- matrix(data = 0, nrow = horizon, ncol = length(irf_list))
    stats_temp <- as.data.frame(matrix(0, nrow = horizon, ncol = 3))
    colnames(stats_temp) <- c("lower", "median", "upper")

    # Get impulse-response for country
    for (draw in seq_along(irf_list)) {
      temp0 <- irf_list[[draw]]
      df_temp[, draw] <- rev(temp0[1:horizon, ncol]) ##/!\
    }

    # Get stats
    stats_temp[, 1] <- apply(df_temp, 1, function(x) quantile(x, probs = 0.16))
    stats_temp[, 2] <- apply(df_temp, 1, median)
    stats_temp[, 3] <- apply(df_temp, 1, function(x) quantile(x, probs = 0.84))

    # Output
    stats_df[[country]] <- stats_temp
  }

  return(stats_df)
}

IRF_processing <- function(A, B, wm) {
  # ==========================================
  # Purpose: computes the irf for a set of countries
  # Parameters:
  #           A: matrices A
  #           B: matrices B
  #           wm: trading weight matrix
  # Returns: for each country a list of irfs
  # ==========================================

  # List of countries
  countries <- as.character(row.names(wm))

  # Initialize output
  irf_total <- list()

  # Initialize counts
  c <- 1
  c2 <- 1

  # Loop on countries
  for (country in countries) {

    A_i <- list()
    B_i <- list()
    irf_country <- list()

    # For each draw, get the country specific matrices A and B
    for (i in seq_along(A)) {
      A_temp <- A[[i]]
      B_temp <- B[[i]]

      A_i[[i]] <- A_temp[c:(c + 2), c:(c + 2)]
      B_i[[i]] <- B_temp[c:(c + 2), c2:(c2 + 7)]

    }

    # Store IRF
    irf_total[[country]] <-  IRF_country(A_i, B_i, 2, 5)

    # Increment count for next country
    c <- c + 3
    c2 <- c2 + 8

  }

  return(irf_total)
}

IRF_country <- function(A_draw, B_draw, nblag, horizon) {
  # ==========================================
  # Purpose: Computes the IRF for each country
  # Parameters:
  #           A_draw: matrices estimated using the MH algo
  #           B_draw: matrices B associated with the drawn matrices A
  #           nblag: number of lags in the model
  #           horizon: time period on which to do the irf
  # Returns: the irf each draw of the country
  # ==========================================

  # Initialization
  irf_results <- list()
  #Two sets of IRF: real GDP growth on precipitation or temperature shocks

  # Loop the draws (compute one irf by draw)
  for (d in seq_along(A_draw)) {

    # Get matrices A and B of the draw
    a <- A_draw[[d]]
    b <- B_draw[[d]]

    nvar <- ncol(a) # number of variables in system
    nrow <- horizon + 1

    # Outpout matrix
    irf_draw <- matrix(data = 0, nrow = (horizon + nblag), ncol = (nvar * nvar))

    # Compute reduced form coefficient: Pi (eq. 11)
    h <- solve(a)
    pi_coeff <- t(h %*% b)
    # reminder: col -> temp.1, prec.1, gdp.1, temp.2, prec.2, gdp.2, intercept, x*

    # Compute irf for each variable
    for (i in seq_len(nvar)) {

      # Initialize last row of IRF
      irf_draw[nrow, (nvar * (i - 1) + 1):(nvar * ((i - 1) + 1))] <- h[i, ]

      # Fill the table backward
      for (j in seq((nrow - 1), 1, -1)) {

        temp0 <- t(irf_draw[j:(j + nblag - 1),(nvar * (i - 1) + 1):(nvar * ((i - 1) + 1))])
        temp1 <- as.matrix(c(temp0[, 1], temp0[, 2]))
        temp2 <- pi_coeff[1:(nvar * nblag), ]

        irf_draw[j, (nvar * (i - 1) + 1):(nvar * ((i - 1) + 1))] <- t(temp1) %*% temp2

      }
    }

    # Add irf to output
    irf_results[[d]] <- irf_draw
  }

  return(irf_results)
}

draw_B <- function(A_revised, X, Y, ar_omega, p) {
  # ==========================================
  # Purpose: computes matrices B associated to the drawn matrices A
  # Parameters:
  #           A_revised: the drawn matrices A
  #           X, Y: exog and endog variables
  #           ar_omega: matrix of variance-covariance from the AR(1) regression
  #           p: number of lags
  # Returns:
  #         B_output: matrices B associated to the draw matrices A
  # ==========================================

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
  # ==========================================
  # Purpose: generate a vector of draws from a normal distribution
  # Parameters:
  #           mean: mean of the normal distribution
  #           cov_var: variance-covariance matrix of the normal distribution
  # Returns:
  #        output: a vector of draws
  # ==========================================

  # Output vector
  output <- c()

  for (i in seq_along(mean)) {

    draw <- rnorm(1, mean = mean[i], sd = sqrt(cov_var[i, i]))

    output <- c(output, draw)
  }

  return(output)
}

draw_D <- function(A_draw, zeta_star_draw, kappa, nb_obs, ar_omega) {
  # ==========================================
  # Purpose: draws all matrices D associated with the matrices A
  # Parameters:
  #           A_draw: matrices A drawed using the MH algo
  #           zeta_star_draw: the zeta_star associated to the A matrices
  #           kappa: constant term
  #           nb_obs: number of observations in the model
  #           ar_omega: matrix of variance-covariance from the AR(1) regression
  # Returns:
  #         output_D: matrices D associated to the draw
  # ==========================================

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
  # ==========================================
  # Purpose: Draws multiple values from an inverse gamma
  #          they all have the same mean but a different variance
  # Parameters:
  #           kappa_star: mean of the inverse gamma
  #           tau_star: vector of variances of the inverse gamma
  # Returns:
  #         draws: a vector of draws from the inverse gamma
  # ==========================================

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
  # ==========================================
  # Purpose: Revise a matrix A with the information set on H
  # Parameters:
  #           A: matrix of structural coefficients to update
  #           mu: weight of matrix H_F
  #           weight_matrix: trading weight matrix used to compute H_F
  # Returns:
  #         A_updated: the new matrix A, updated with infos from H
  # ==========================================

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

metropolis_hastings_algorithm <- function(iter, burnin,
                                          A_ini, X, Y, kappa, omega, ar_omega,
                                          prior_spec, weight_matrix, pU) {
  # ==========================================
  # Purpose: lauches the MH algorithm to draw a sample of matrices A
  # Parameters:
  #           iter, burnin: number of iterations and burnin
  #           A_ini: initial matrix A (source: GVARX)
  #           X, Y: exog and endog variables
  #           kappa: constant term
  #           omega: variance-covariance matrix of the model reduced form
  #           ar_omega: variance-covariance matrix of AR(1) regression
  #           prior_spec: information on the prior of some coefficients of A
  #           weight_matrix: trading weight matrix
  #           pU: coefficients of A that have to be drawn
  # Returns:
  #         output_A: list storing the draws of A
  #         output_zeta: list storing the associated zetas
  # ==========================================

  # Initialization - GVARX is the starting point
  A_old <- A_ini
  regr <- compute_tilded_regr(A_old, X, Y, ar_omega, p = 2)
  zeta_old <- regr[[1]]
  A_old_post <- posterior_A(A_old, X, Y, kappa, zeta_old, omega, ar_omega,
                            prior_spec, weight_matrix)

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
    A_new_post <- posterior_A(A_new, X, Y, kappa, zeta_new, omega, ar_omega,
                              prior_spec, weight_matrix)

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
  # ==========================================
  # Purpose: draws a matrix A after a normal shock
  # Parameters:
  #           A: matrix A, starting point for the draw
  #           pU: coefficients of matrix A to update
  #           scale: scale of the shock
  #           wm: trading weight matrix
  # Returns:
  #           A_updated: the new A matrix
  # ==========================================

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

posterior_A <- function(A, X, Y, kappa, zeta_star, 
                        omega, ar_omega, prior_spec, weight_matrix) {
  # ==========================================
  # Purpose: function that computes the posterior of A (log)
  # Parameters:
  #           A, X, Y: matrix of structural coeff, exog and endog variables
  #           kappa: constant term
  #           zeta_star: sum of squared residuals of tilded regression
  #           omega: variance-covariance matrix of the model reduced form
  #           ar_omega: variance-covariance matrix of the AR(1) regression
  #           prior_spec: parameters of the student distribution for some coefficients of A
  #           weight_matrix: matrix containing trading weights
  # Returns:
  #         posterior: the posterior of A
  # ==========================================

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
  # ==========================================
  # Purpose: computes the log-likelihood associated with a matrix A
  # Parameters:
  #           A: structural coefficients
  #           X: exogeneous variables
  #           Y: endogeneous variables
  #           kappa: constant term
  #           zeta_star: sum of squared residuals of tilded regression
  #           omega: variance-covariance matrix of the model reduced form
  #           ar_omega: variance-covariance matrix of the AR(1) regression
  # Returns:
  #         like: the log-likelihood of matrix A
  # ==========================================

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
  # ==========================================
  # Purpose: computes the tilded regression detailed on page 24
  #          of the article (regression with augmented vectors)
  # Parameters:
  #           A: matrix of structural coefficients
  #           X: matrix of exogeneous variables
  #           Y: matrix of endogeneous variables
  #           p: number of lag
  # Returns:
  #         zeta: sum of squared residuals of the regression
  #         m_star: mean
  #         M_star: matrix of variance-covariance
  # ==========================================

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
  # ==========================================
  # Purpose: compute the log of the prior of matrix A
  # Parameters:
  #           A: matrix A with the coefficients of the VAR regression
  #           prior_spec: parameters of the distribution for the
  #                       coefficients a_gdp_temp and a_gdp_prec
  #           weight_matrix: matrix containing trading weights to compute H_F
  # Returns:
  #        prior: the prior of matrix A (still a log value),
  #               which is the sum of each country's priors
  # ==========================================

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
  # ==========================================
  # Purpose: draws a value from the density of a student distribution
  # Parameters:
  #           x: value of parameter
  #           mu: mode
  #           nu: degree of freedom
  #           sigma: scale
  # Returns:
  #       prior: prior of the parameter inputed
  # ==========================================

  # Load library
  library(LaplacesDemon)

  # Compute prior
  prior <- dst(x = x, mu = mu, sigma = sigma, nu = nu, log = FALSE)

  return(prior)
}

aromega_computation <- function(y, p) {
  # ==========================================
  # Purpose: Computes the variance-covariance matrix of an Ar(1) regression
  # Parameters:
  #           Y: vector containing the endogeneous variables
  #           p: number of lags
  # Returns:
  #        var_cov: the variance-covariance matrix
  # ==========================================


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