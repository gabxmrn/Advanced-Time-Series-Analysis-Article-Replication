source("code/models/linear_regr.R")

compute_ic <- function(data, country, max_lag) {
  # ==========================================
  # Purpose: compute the AIC, BIC and HQN for a VAR with constant
  # Parameters:
  #           data: dataframe containing data (columns = dates)
  #           country: country of the dataframe
  #           max_lag: maximum lag to test
  # Returns: nothing
  # ==========================================

  # Data preparation
  data <- t(data)
  k <- ncol(data) # Number of variables
  n <- nrow(data) # Number of observations

  # Empty dataframe for result
  results <- data.frame(matrix(nrow = max_lag, ncol = 3))
  colnames(results) <- c("AIC", "BIC", "HQN")

  # Loop on the number of lags
  for (p in 1:max_lag) {

    # Create lag data
    y <- data[(p + 1):n, ]
    x <- do.call(cbind, lapply(1:p, function(lag) data[(p + 1 - lag):(n - lag), ]))

    # Compute linear regression
    regr <- linear_reg(y, x, constant = TRUE)
    res <- regr$residuals # Residuals of model
    n_res <- nrow(res)
    sigma_res <- crossprod(res) / n_res # Covariance matrix of residuals

    # # Compute log-likelihood
    # loglikelihood <- - 0.5 * n_res * k * log(2 * pi) -
    #   0.5 * n_res * log(det(sigma_res)) -
    #   0.5 * sum(diag(res %*% solve(sigma_res) %*% t(res)))

    # # Compute information criterion
    # # Result: all countries should have 5 lags according to the 3 criterion
    # results[p, 1] <- 2 * (k * p + 1) - 2 * loglikelihood
    # results[p, 2] <- (k * p + 1) * log(n_res) - 2 * loglikelihood
    # results[p, 3] <- 2 * (k * p + 1) * log(log(n_res)) - 2 * loglikelihood

    # Other method -> uses a proxy for log-likelihood (see vars package)
    # Does not give same results -> closer to the ones in the article (majority of countries with lag = 1)
    results[p, 1] <- log(det(sigma_res)) + (2 / n_res) * (p * k^2 + k)
    results[p, 2] <- log(det(sigma_res)) + (log(n_res) / n_res) * (p * k^2 + k)
    results[p, 3] <- log(det(sigma_res)) + (2 * log(log(n_res)) / n_res) *
      (p * k^2 + k)
  }

  # Display results
  cat("Information Criterion - Result for ", country, "\n")
  print(results)

  # Print smallest lag for each information criteria
  cat("Smallest lag - AIC: ", which.min(results[, 1]), "\n")
  cat("Smallest lag - BIC: ", which.min(results[, 2]), "\n")
  cat("Smallest lag - HQN: ", which.min(results[, 3]), "\n")

}