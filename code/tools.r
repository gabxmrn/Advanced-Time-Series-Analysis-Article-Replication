data_importation <- function(path, sheet_name) {
  # ==========================================
  # Function: data_importation
  # Purpose: Import data from an Excel file and format it as a dataframe
  #          with rows being countries and columns years.
  # Parameters:
  #            path (string): path to the excel file.
  #            sheet_name (string): name of the sheet in the file.
  # Returns:
  #            df (dataframe): a dataframe with rows being countries
  #                            and columns years.
  # ==========================================
  library(readxl)
  df <- read_excel(path, sheet = sheet_name,
                   col_names = TRUE) # Import data from excel
  df <- as.data.frame(df) # Convert type to dataframe
  rownames(df) <- df[[1]] # Set index to country names
  df <- df[-1] # Delete column with country names
  return(df)
}

climate_variable_treatment <- function(df) {
  # ==========================================
  # Function: climate_variable_treatment
  # Purpose: addjusts climate variables by subtracting a 30-year moving
  #          average from the value of each year.
  # Parameters:
  #            df: input dataframe where rows are countries
  #                and columns are years.
  # Returns: new dataframe containing adjusted values.
  # ==========================================
  df_result <- data.frame(row.names = rownames(df))

  years <- as.numeric(colnames(df)) # Get possible years in the dataframe

  # Iterate over each year
  for (year in years) {
    # Compute start year for the moving average
    start_year <- year - 31

    if (start_year %in% years) {
      # Average on the last 30 years (ending at t-1)
      mov_av <- rowMeans(df[, which(years >= start_year & years < year)],
                         na.rm = TRUE)
      # Substracts the moving average to the value of year t
      df_result[[as.character(year)]] <- df[, as.character(year)] - mov_av
    }
  }

  return(df_result)
}