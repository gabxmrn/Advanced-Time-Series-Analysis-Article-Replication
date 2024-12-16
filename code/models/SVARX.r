source("code/utils/information_criteria.R")

svarx_main <- function(temperatures, precipitations, gdp, foreign_variable) {
  # ==========================================
  # Purpose:
  # Parameters:
  # Returns:
  # ==========================================

  countries <- as.character(row.names(temperatures))

  for (country in countries) {

    # Creation of matrix y
    y <- rbind(
      temperatures[rownames(temperatures) == country, ],
      precipitations[rownames(precipitations) == country, ],
      gdp[rownames(gdp) == country, ]
    )
    rownames(y) <- c("temperatures", "precipitations", "gdp")

    # Lag selection - information criterion
    if (country == "Austria") {
      compute_ic(y, country, 5)
    }
    # compute_ic(y, country, 5) # compute ic

    # Foreign variable for country i
    x_star <- foreign_variable[rownames(foreign_variable) == country]

    # Creation of matrix x (foreign variable matrix)
    nblag <- 2 # Forced value from article
    x <- svarx_fvm(y, x_star, p = nblag)

    # Dropped first values of y due to lag
    n <- ncol(y)
    y <- y [, (nblag + 1):n]

    # SVAR estimation
    if (country == "Austria") {
      svarx_estm(y, x, p = nblag)
    }
  }

}

svarx_estm <- function(y, x, p) {

}

svarx_fvm <- function(y, x_star, p) {
  # ==========================================
  # Purpose: created the vector x for the nested SVARX composed of
  #          the lags of the country, an intercept and foreign variable
  # Parameters:
  #           y: endogenous data set
  #           weight_matrix: weight matrix for the foreign variable
  #           p: number of lags for the SVARX
  # Returns:
  #        x: exogenous variables for the SVARX model
  # ==========================================

  # Add lagged y
  y <- t(y)
  n <- nrow(y)
  x <- do.call(cbind, lapply(1:p, function(lag) y[(p + 1 - lag):(n - lag), ]))
  x <- t(x)
  lag <- 1
  for (i in 1:p) {
    rownames(x)[lag] <- paste0("temp.", i) # Rename row
    rownames(x)[lag + 1] <- paste0("prec.", i) # Rename row
    rownames(x)[lag + 2] <- paste0("gdp.", i) # Rename row
    lag <- lag + 3
  }

  # Add intercept and foreign variable
  x <- rbind(x, 1, x_star[2:(n - 1)])
  rownames(x)[p * 3 + 1] <- "intercept" # Rename row
  rownames(x)[nrow(x)] <- "x.star" # Rename row

  # Reindex dataframe
  colnames(x) <- NULL

  return(x)
}