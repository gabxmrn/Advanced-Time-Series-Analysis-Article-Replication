source("code/data_management.R")
source("code/stationnarity.R")

########## Data Importation ##########

# Import data from Excel sheets
temperatures_i <- data_importation("code/data/data.xlsx", "temperatures")
precipitations_i <- data_importation("code/data/data.xlsx", "precipitations")
gdp_i <- data_importation("code/data/data.xlsx", "GDP")
weight_maxtrix <- data_importation("code/data/trade balance.xlsx", "Feuil1")

########## Data Treatment ##########

# Climate variables: substract a 30Y moving average from the value of each year
temperatures <- climate_variable_treatment(temperatures_i)
precipitations <- climate_variable_treatment(precipitations_i)

# GDP : log difference of the original series
gdp <- t(apply(log(gdp_i), 1, diff))
gdp <- gdp[, as.character(1991:2019)]

# Foreign exogenous variable: product of weight matrix and GDP growth variable

########## Plotting series ##########
# plot_series(temperatures, "temperatures", "done")
# plot_series(precipitations, "precipitations", "done")
# plot_series(gdp, "gdp", "done")

########## Testing stationnarity ##########
kpss_temp <- kpss(temperatures, 0.01, "level")
print(kpss_temp)

# Rq: les séries temperatures/precipitations/gdp sont statio
#     au seuil de 1% (avec constante, sans tendance)