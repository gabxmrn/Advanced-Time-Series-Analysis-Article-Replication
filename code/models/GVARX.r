source("code/utils/information_criteria.R")

gvarx_main <- function(temperatures, precipitations, gdp, trade_balance) {
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

    # Lag selection -> force_lag or call function
    compute_ic(y, country, 5) # compute ic
    nb_lags <- 2 # force lag number to article one

    # Creation of matrix x (foreign variable matrix)
    # x <- gvarx_ft(y, gdp, trade_balance)

  }



}

gvarx_ft <- function(y, gdp, trade_balance, nb_lags) {

}