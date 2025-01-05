


bs_main <- function(a, b, y, x, u, prior_specifications) {

  prior_a(prior_specifications)

}


posterior_a <- function() {

}

prior_a <- function(prior_spec) {

  library(LaplacesDemon)

  countries <- as.character(row.names(prior_spec))

  # Loop on countries
  for (country in countries) {

    # Get model parameters - a_gdp_temp and a_gdp_prec
    mode_gt_gp <- prior_spec[country, 1]
    scale_gt_gp <- prior_spec[country, 2]
    dof_gt_gp <- prior_spec[country, 3]

    # Get model parameters - a_temp_gdp and a_prec_gdp
    mode_tg_pg <- 0
    scale_tg_pg <- 0.01
    dof_tg_pg <- 3

    # Priors of the elements of A
    prior_a_t_g <- dst(x = seq(-0.06, 0.06, length.out = 10000), mu = mode_tg_pg, sigma = scale_tg_pg, nu = dof_tg_pg, log = FALSE)
    prior_a_p_g <- dst(x = seq(-0.06, 0.06, length.out = 10000), mu = mode_tg_pg, sigma = scale_tg_pg, nu = dof_tg_pg, log = FALSE)
    prior_a_g_t <- dst(x = seq(-0.06, 0.06, length.out = 10000), mu = mode_gt_gp, sigma = scale_gt_gp, nu = dof_gt_gp, log = FALSE)
    prior_a_g_p <- dst(x = seq(-0.06, 0.06, length.out = 10000), mu = mode_gt_gp, sigma = scale_gt_gp, nu = dof_gt_gp, log = FALSE)

    # Prior of H_D

    # Prior of H_F

  }

  plot(prior_a_g_p)

  # Product of all p(a_i)

  # Return prior of A

}