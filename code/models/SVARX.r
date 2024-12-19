source("code/utils/information_criteria.R")
source("code/models/linear_regr.R")

svarx_main <- function(temperatures, precipitations, gdp, foreign_variable) {
  # ==========================================
  # Purpose: Estimate country-specific VAR and stack resulting matrix
  # Parameters: 4 dataframes containing model datas
  # Returns:
  #         A: global structural matrix
  #         B: global coefficient matrix
  #         Y: global Y matrix
  #         X: global  X matrix
  #         U: global residual matrix
  # ==========================================

  # Load libraries
  library("Matrix")

  countries <- as.character(row.names(temperatures))
  n <- ncol(temperatures)
  nblag <- 2 # Force lag value

  # Output creation
  a <- matrix(nrow = 0, ncol = 0)
  b <- matrix(nrow = 0, ncol = 0)
  u <- matrix(nrow = 0, ncol = (n - nblag))
  y <- matrix(nrow = 0, ncol = (n - nblag))
  x <- matrix(nrow = 0, ncol = (n - nblag))

  for (country in countries) {

    # Creation of matrix y
    y_i <- rbind(
      temperatures[rownames(temperatures) == country, ],
      precipitations[rownames(precipitations) == country, ],
      gdp[rownames(gdp) == country, ]
    )
    rownames(y_i) <- c("temperatures", "precipitations", "gdp")

    # Lag selection - information criterion
    # compute_ic(y, country, 5) # compute ic

    # Foreign variable for country i
    x_star <- foreign_variable[rownames(foreign_variable) == country]

    # Creation of matrix x (foreign variable matrix)
    x_i <- svarx_fvm(y_i, x_star, p = nblag)

    # Dropped first values of y due to lag
    y_i <- y_i[, (nblag + 1):n]

    # SVAR estimation
    temp <- svarx_estm(y_i, x_i, p = nblag)
    a_i <- temp[[1]]
    b_i <- temp[[2]]
    u_i <- temp[[3]]

    # Stack matrix
    a <- bdiag(a, a_i)
    b <- bdiag(b, b_i)

    # Output - Matrix Y
    y <- rbind(y, y_i)
    rownames(y)[nrow(y) - 2] <- paste0(country, ".temp")
    rownames(y)[nrow(y) - 1] <- paste0(country, ".prec")
    rownames(y)[nrow(y)] <- paste0(country, ".gdp")

    # Output - Matrix X
    x <- rbind(x, x_i)
    count <- 2
    for (i in nblag:1) {
      rownames(x)[nrow(x) - count - 2] <- paste0(country, ".temp.", i)
      rownames(x)[nrow(x) - count - 1] <- paste0(country, ".prec.", i)
      rownames(x)[nrow(x) - count] <- paste0(country, ".gdp.", i)
      count <- count + 3
    }
    rownames(x)[nrow(x) - 1] <- paste0(country, ".intercept")
    rownames(x)[nrow(x)] <- paste0(country, ".x_star")

    # Output - Matrix U
    u <- rbind(u, u_i)
    rownames(u)[nrow(u) - 2] <- paste0(country, ".temp")
    rownames(u)[nrow(u) - 1] <- paste0(country, ".prec")
    rownames(u)[nrow(u)] <- paste0(country, ".gdp")

  }

  colnames(x) <- colnames(y) # Rename columns of X

  # Output - Matrix A and B
  a_df <- as.data.frame(as.matrix(a))
  b_df <- as.data.frame(as.matrix(b))

  count2 <- 1
  count3 <- 1
  for (country in countries) {
    rownames(a_df)[count2] <- paste0(country, ".temp")
    rownames(a_df)[count2 + 1] <- paste0(country, ".prec")
    rownames(a_df)[count2 + 2] <- paste0(country, ".gdp")

    colnames(a_df)[count2] <- paste0(country, ".temp")
    colnames(a_df)[count2 + 1] <- paste0(country, ".prec")
    colnames(a_df)[count2 + 2] <- paste0(country, ".gdp")

    rownames(b_df)[count2] <- paste0(country, ".temp")
    rownames(b_df)[count2 + 1] <- paste0(country, ".prec")
    rownames(b_df)[count2 + 2] <- paste0(country, ".gdp")

    for (count4 in 1:nblag) {
      colnames(b_df)[count3] <- paste0(country, ".temp.", count4)
      colnames(b_df)[count3 + 1] <- paste0(country, ".prec.", count4)
      colnames(b_df)[count3 + 2] <- paste0(country, ".gdp.", count4)
      count3 <- count3 + 3
    }
    colnames(b_df)[count3] <- paste0(country, ".intercept")
    colnames(b_df)[count3 + 1] <- paste0(country, ".x_star")

    count3 <- count3 + 2
    count2 <- count2 + 3
  }

  # Return matrix A, B, Y, X and U
  return(list(a_df, b_df, y, x, u))
}

svarx_estm <- function(y, x, p) {
  # ==========================================
  # Purpose: Estimates the SVARX model for country i
  # Parameters:
  #           y: endogenous variables
  #           x: exogenous variables
  #           p: number of lags in the model
  # Returns:
  #         A: structural matrix for country i
  #         B: coefficient matrix for country i
  #         U: residual matrix for country i
  # ==========================================

  # Model parameters
  k <- nrow(y) # Number of exogenous variables
  n <- ncol(y) # Number of observations

  # Output variable
  a <- diag(k)
  rownames(a) <- cbind("temp", "prec", "gdp")
  colnames(a) <- cbind("temp", "prec", "gdp")

  b <- matrix(nrow = k, ncol = (p * k + 2))
  rownames(b) <- cbind("temp", "prec", "gdp")
  colnames(b) <- rownames(x)

  u <- matrix(nrow = k, ncol = n)
  rownames(u) <- cbind("temp", "prec", "gdp")
  colnames(u) <- colnames(y)

  # VAR estimation without constant (already included in x)
  for (i in 1:k) {

    colnames(x) <- colnames(y)

    # Adding data for matrix A parameters computation
    if (i == 3) {
      x_a <- rbind(y[1, ], y[2, ], x)
      rownames(x_a)[1] <- "temp"
      rownames(x_a)[2] <- "prec"
    } else {
      x_a <- rbind(y[3, ], x)
      rownames(x_a)[1] <- "gdp"
    }

    # Linear regression
    regr <- linear_reg(t(y[i, ]), t(x_a), constant = FALSE)

    # Output - variable A
    if (i == 1) {
      a[1, 3] <- - regr$coefficients[1]
    } else if (i == 2) {
      a[2, 3] <- - regr$coefficients[1]
    } else {
      a[3, 1] <- - regr$coefficients[1]
      a[3, 2] <- - regr$coefficients[2]
    }

    # Output - variable B
    n_temp <- nrow(regr$coefficients)
    if (i == 3) {
      b[i, ] <- t(regr$coefficients[3:n_temp, ])
    } else {
      b[i, ] <- t(regr$coefficients[2:n_temp, ])
    }

    # Output - variable u
    u[i, ] <- t(regr$residuals)
  }

  return(list(a, b, u))
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