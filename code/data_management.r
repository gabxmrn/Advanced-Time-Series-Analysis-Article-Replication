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

plot_series <- function(df, t_series, process) {
  # ==========================================
  # Function: plot_series
  # Purpose: Display a time series (temperatures, precipitations
  #          or gdp) for each country on a 5x5 grid
  # Parameters:
  #           df: dataframe containing the data
  #           t_series: name of the time series
  #           process: if it's the original series or if it has been transformed
  # Returns:
  #           nothing, plots the 5x5 grid of graphs
  # ==========================================}

  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(reshape2)

  # Check params
  if (tolower(process) != "none" && tolower(process) != "done") {
    stop("Process should be either none or done.")
  }

  # Y-axis title
  if (tolower(t_series) == "temperatures") {
    y_axis_legend <- "Temperatures (°C)"
  } else if (tolower(t_series) == "precipitations") {
    y_axis_legend <- "Precipitations (mm)"
  } else if (tolower(t_series) == "gdp") {
    y_axis_legend <- "GDP (mil. 2017US$)"
  } else {
    stop("Variable type should be temperatures, precipitations or gdp.")
  }

  # Data management before plotting
  df_t <- as.data.frame(t(df))
  colnames(df_t) <- rownames(df)
  df_t$Year <- as.numeric(row.names(df_t))

  df_t_l <- melt(df_t, id.vars = "Year", variable.name = "Country",
                 value.name = "Value")

  # Plotting the graph for each country
  plots <- lapply(unique(df_t_l$Country), function(country) {
    ggplot(df_t_l[df_t_l$Country == country, ], aes(x = Year, y = Value)) + # nolint
      geom_line(color = "navyblue") +
      labs(title = country, x = "Year", y = y_axis_legend) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)
      )
  })

  # Visualisation of all graphs
  if (tolower(t_series) == "temperatures") {
    if (tolower(process) == "none") {
      grid_title <- textGrob("Average Yearly Temperature",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    } else {
      grid_title <- textGrob("Average Yearly Temperature post-processing",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    }
  } else if (tolower(t_series) == "precipitations") {
    if (tolower(process) == "none") {
      grid_title <- textGrob("Average Yearly Precipitations",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    } else {
      grid_title <- textGrob("Average Yearly Precipitations post-processing",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    }
  } else if (tolower(t_series) == "gdp") {
    if (tolower(process) == "none") {
      grid_title <- textGrob("Real GDP at constant 2017 national prices",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    } else {
      grid_title <- textGrob("Real GDP at constant 2017 USD post-processing",
                             gp = gpar(fontsize = 16, fontface = "bold"))
    }
  }
  grid.arrange(grobs = plots, ncol = 5, top =  grid_title)
}