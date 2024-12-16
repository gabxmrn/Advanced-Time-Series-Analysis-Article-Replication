source("code/models/linear_regr.R")

kpss <- function(df, alpha, level_or_trend) {
  # ==========================================
  # Purpose: Tests stationnarity for a dataframe
  #          Null hypothesis: data is stationnary
  # Parameters:
  #           df: dataframe to test, countries must be rows
  #           alpha: level of significance for results display
  #           level_or_trend: test if data is stationnary around
  #                           a constant or a constant + trend
  # Returns:
  #           A dataframe with the p-values and the conclusion of the test
  # ==========================================

  # Check alpha value
  if (alpha >= 1) {
    stop("Alpha should be between 0 and 1.")
  }

  # Check level or trend parameter
  if (level_or_trend != "trend" && level_or_trend != "level") {
    stop("Series should be stationnary around level or trend. 
          Please enter level or trend.")
  }

  # Creation of empty dataframe for results
  df_result <- data.frame(row.names = rownames(df))
  df_result[, "p-values"] <- NA
  df_result[, "stationnary"] <- NA

  countries <- as.character(rownames(df))

  # KPSS Test for each country
  for (country in countries) {

    # Select data for the country
    temp_df <- t(df[country, ])
    n <- length(temp_df)

    # Table of probabilities associated to critical values
    tablep <- c(0.01, 0.025, 0.05, 0.1)

    # Linear regression of the model
    if (level_or_trend == "trend") {
      t <- 1:n
      regr <- linear_reg(y = temp_df, x = t) # linear regression with trend
      stat_table <- c(0.216, 0.176, 0.146, 0.119) # Critical values
    } else {
      regr <- linear_reg(y = temp_df) # Linear regression with constant
      stat_table <- c(0.739, 0.574, 0.463, 0.347) # Critical values
    }

    res <- regr$residuals # Residuals of model

    # Compute test statistic
    s <- cumsum(res) # Cumulative sum of residuals
    eta <- sum(s^2) / (n^2) # Partial sum
    s2 <- sum(res^2) / n # Average square residuals

    lm_stat <- eta / s2  # Test statistic

    # Compute p-value
    p_val <- approx(stat_table, tablep, lm_stat, rule = 2)$y

    # Fill dataframe with pvalues and test conclusion
    df_result[country, 1] <- round(p_val, 2)

    if (p_val < alpha) {
      df_result[country, 2] <- "No"
    } else {
      df_result[country, 2] <- "Yes"
    }
  }

  return(df_result)
}