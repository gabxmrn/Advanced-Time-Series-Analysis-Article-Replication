source("code/utils/data_importation.R")
source("code/utils/stationnarity.R")
source("code/models/SVARX.R")
source("code/models/BSVAR.R")

########## Data Importation ##########

temperatures_i <- data_importation("code/data/data.xlsx", "temperatures")
precipitations_i <- data_importation("code/data/data.xlsx", "precipitations")
gdp_i <- data_importation("code/data/data.xlsx", "GDP")
weight_matrix <- data_importation("code/data/trade balance.xlsx", "Feuil1")

########## Data Treatment ##########

# Climate variables: substract a 30Y moving average from the value of each year
temperatures <- climate_variable_treatment(temperatures_i)
precipitations <- climate_variable_treatment(precipitations_i)

# GDP : log difference of the original series
gdp <- t(apply(log(gdp_i), 1, diff))
gdp <- as.data.frame(gdp[, as.character(1991:2019)])

# Foreign exogenous variable computation
foreign_var <- as.matrix(weight_matrix) %*% as.matrix(gdp)

# Import data for prior
prior_spec <- data_importation("code/data/priors_specifications.xlsx", "Feuil1")

########## Plotting series ##########
# plot_series(temperatures, "temperatures", "done")
# plot_series(precipitations, "precipitations", "done")
# plot_series(gdp, "gdp", "done")

########## Testing stationnarity ##########
# kpss_temp <- kpss(precipitations, 0.01, "level")
# print(kpss_temp)

########## Model ##########

svar_results <- svarx_main(temperatures, precipitations, gdp, foreign_var)
bsvarx_results <- bs_main(svar_results, prior_spec, weight_matrix)