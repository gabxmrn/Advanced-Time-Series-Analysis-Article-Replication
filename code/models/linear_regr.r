linear_reg <- function(y, x =  NULL, constant = TRUE) {
  # ==========================================
  # Purpose: computes the linear regression for the models y ~ 1
  #          and y ~ t with t being a trend
  # Parameters:
  #           y: endogenous variable
  #           x: constant or constant + trend
  # Returns:
  #           beta: beta of the model
  #           fitted_values: fitted values of the model
  #           residuals: residuals of the model
  # ==========================================

  if (is.null(x)) {
    # Case of an intercept only model
    x <- matrix(rep(1, length(y)), ncol = 1)
  } else {
    x <- as.matrix(x)
    # Add intercept
    if (constant) {
      x <- cbind(1, x)
    }
  }

  # Beta: (t(x)x)^(-1)t(x)y
  beta <- solve(t(x) %*% x) %*% t(x) %*% y

  # Fitted values
  fitted_values <- x %*% beta

  # Residuals
  residuals <- y - fitted_values

  # Return results
  return(list(
    coefficients = beta,
    fitted_values = fitted_values,
    residuals = residuals
  ))

}
