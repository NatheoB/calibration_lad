initialise_models <- function(exp_design,
                              data_sensors,
                              data_plots,
                              data_calib,
                              data_rad,
                              output_lad_method1,
                              species2calib) {
  
  
  ## Create the Bayesian setup for each model ----
  model_setups <- vector("list", nrow(exp_design))
  for (i in 1:nrow(exp_design)) {
    
    model_setups[[i]] <- initialise_model(exp_design[i,],
                                          data_sensors, 
                                          data_plots,
                                          data_calib, 
                                          data_rad,
                                          output_lad_method1,
                                          species2calib)
  }
  
  
  ## Return the model setups ----
  return(list(
    "setups" = model_setups,
    "exp_design" = exp_design
  ))
}



initialise_model <- function(mod_design,
                             data_sensors, 
                             data_plots,
                             data_calib, 
                             data_rad,
                             output_lad_method1,
                             species2calib) {
  
  # BayesianTools consider global variables during the MCMC sampling
  # Thus, we need to create an environment with all the functions and variables needed for MCMC sampling
  
  
  
  # SUB-FUNCTIONS ----
  
  ## FUNCTION TO RUN SAMSARALIGHT ----
  run_sl_standXlad <- function(data_plots,
                               data_calib, 
                               data_rad,
                               site, 
                               lad_sites,
                               lad_species,
                               lad_dbh,
                               lad_compet) {
    
    # Create table of intercepts for each species
    # Here the order of the species are really important
    # It is the same as in the species2calib and in p_lad_intercepts
    sp_intercepts <- data.frame(
      species_calib = species2calib,
      intercept = lad_species,
      dbh_coef = lad_dbh,
      compet_coef = lad_compet
    )
    
    # Get the random effect of the site
    # Here, the order is really important, and is the same between site_names and lad_sites vectors
    site_intercept <- lad_sites[which(site_names == site)]
    
    # Get tree dataset and set the LAD value for all trees
    # Here, only add the site random effect, as the origin effect in already included in the estimation of the site effect (hierarchical structure)
    tmp_trees <- data_calib[[site]]$trees %>% 
      dplyr::left_join(sp_intercepts, by = "species_calib") %>% 
      dplyr::mutate( 
        crown_lad = site_intercept + intercept + dbh_cm * dbh_coef + balh_m2ha * compet_coef
      )
    
    # Get the plot info
    tmp_plot <- data_plots %>% 
      dplyr::filter(name == site)
    
    # Run SamsaraLight
    tmp_out_samsalight <- 
      sl_run(tmp_trees, 
             data_rad[[site]],
             sensors = data_calib[[site]]$sensors, 
             sensors_only = TRUE,
             latitude = tmp_plot$latitude, 
             slope = tmp_plot$slope, 
             aspect = tmp_plot$aspect, 
             north_to_x_cw = tmp_plot$northToX,
             start_day = 121, 
             end_day = 273,
             cell_size = data_calib[[site]]$info$cell_size, 
             n_cells_x = data_calib[[site]]$info$n_cells_x, 
             n_cells_y = data_calib[[site]]$info$n_cells_y,
             turbid_medium = TRUE,
             trunk_interception = TRUE,
             soc = TRUE,
             height_anglemin = 15,
             direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
             direct_anglestep = 5,
             diffuse_anglestep = 15,
             detailed_output = FALSE)
    
    # Sl sensors output for the plot
    tmp_out_samsalight$output$sensors
    
  }
  
  ## FUNCTION TO RUN SAMSARALIGHT ON ALL STANDS ----
  compute_pacl_residuals <- function(data_sensors,
                                     data_plots,
                                     data_calib,
                                     data_rad,
                                     lad_sites,
                                     lad_species,
                                     lad_dbh,
                                     lad_compet,
                                     print.pb) {
    
    # Output list
    out_residuals_list <- vector("list", length(site_names))
    
    # Initialize progress bar
    i <- 1
    
    if (print.pb) {
      pb <- txtProgressBar(min = 0, max = length(site_names),
                           style = 3, width = 50, char = "=")
    }
    
    # Run for each simulation of the experimental design
    for (site in site_names) {
      
      # Sl sensors output for a given LAD model in a given plot
      tmp_out_sl <- run_sl_standXlad(data_plots,
                                     data_calib, 
                                     data_rad,
                                     site, 
                                     lad_sites,
                                     lad_species,
                                     lad_dbh,
                                     lad_compet)
      
      # Compute mean residuals between all sensors
      out_residuals_list[[i]] <- dplyr::left_join(
        
        # Estimated pacl from virtual sensors
        tmp_out_sl %>%
          dplyr::select(id_sensor, pacl_sl = pacl),
        
        # Measured pacl from field sensor
        data_sensors[[site]] %>%
          dplyr::select(id_sensor = id, pacl_field = PACLtotal),
        
        by = "id_sensor"
      ) %>% 
        
        dplyr::mutate(residuals = pacl_sl - pacl_field) %>% 
        dplyr::select(id_sensor, pacl_sl, pacl_field, residuals)
      
      
      # Update the progress bar
      if (print.pb) setTxtProgressBar(pb, i)
      i <- i+1
    }
    if (print.pb) close(pb)
    
    dplyr::bind_rows(out_residuals_list)
  }
  
  
  
  ## FUNCTION TO CREATE THE SUBVECTOR OF PARAMETERS ----
  split_parameters_vector <- function(p) {
    
    # 1. Base model parameters
    p_sigma <- p[1] # Standard deviation of residuals
    p_sigma_origins <- p[2] # SD of the origin effect
    p_sigma_sites <- p[3] # SD of the site effect
    p_lad_origins <- p[(1:n_origins) + 3] # Origins intercepts
    p_lad_sites <- p[(1:n_sites) + 3 + n_origins] # Site intercepts
    
    
    #### ICI: voir pour moduler selon si random effects !
    ### TIPS: mettre une variable start_species, start dbh...
    
    # 2. Parameters depending on the model definition
    if (mod_design$species_specific) {
      # species-specific intercept
      p_lad_species <- p[(1:n_species) + 3 + n_origins + n_sites]
      
    } else {
      # a single intercept for all species
      p_lad_species <- rep(p[1 + 3 + n_origins + n_sites], n_species)
    }
    
    
    if (mod_design$dbh_effect) {
      
      if (mod_design$species_specific) {
        # dbh effect for each species
        p_lad_dbh <- p[(1:n_species) + (3 + n_origins + n_sites + n_species)]
        
      } else {
        # a single dbh effect for all species
        p_lad_dbh <- rep(p[1 + (3 + n_origins + n_sites + 1)], n_species)
      }
      
    } else {
      # Otherwise, null dbh effect
      p_lad_dbh <- rep(0, n_species)
    } 
    
    
    if (mod_design$compet_effect) {
      
      if (mod_design$species_specific) {
        # compet effect for each species
        # be careful, position on the vector of parameters depends on if dbh effect is considered
        if (mod_design$dbh_effect) {
          p_lad_compet <- p[(1:n_species) + (3 + n_origins + n_sites + 2 * n_species)]
        } else {
          p_lad_compet <- p[(1:n_species) + (3 + n_origins + n_sites + n_species)]
        }
        
      } else {
        # a single compet effect for all species
        # be careful, position on the vector of parameters depends on if dbh effect is considered
        if (mod_design$dbh_effect) {
          p_lad_compet <- rep(p[1 + (3 + n_origins + n_sites + 2)], n_species)
        } else {
          p_lad_compet <- rep(p[1 + (3 + n_origins + n_sites + 1)], n_species)
        }
        
      }
      
    } else {
      # Otherwise, null compet effect
      p_lad_compet <- rep(0, n_species)
    } 
    
    # List with all different group of parameters 
    return(list(
      "sigma" = p_sigma,
      "sigma_origins" = p_sigma_origins,
      "sigma_sites" = p_sigma_sites,
      "origins" = p_lad_origins,
      "sites" = p_lad_sites,
      "species" = p_lad_species,
      "dbh" = p_lad_dbh,
      "compet" = p_lad_compet
    ))
  }
  
  
  ## FUNCTION TO COMPUTE LOG POSTERIOR FROM LIKELIHOODS ----
  compute_log_posterior <- function(p, pointwise = FALSE, print.pb = FALSE) {
    
    # 1. Parameters for the calibration run
    p_list <- split_parameters_vector(p)
    
    
    # 2. Compute the residuals of the given lad parameters
    out_sl <- compute_pacl_residuals(data_sensors,
                                     data_plots,
                                     data_calib, 
                                     data_rad, 
                                     p_list$sites,
                                     p_list$species,
                                     p_list$dbh,
                                     p_list$compet,
                                     print.pb)
    residuals <- out_sl$residuals
    
    
    ## ICI : Compute log likelihood seleuement si random
    
    # 3. Compute the log-likelihood of the origin random effect
    log_likelihood_origins <- dnorm(p_list$origins, mean = 0, sd = p_list$sigma_origins, log = TRUE)
    
    
    
    # 4. Compute the log-likelihood of the hirerarchcal random site effect
    
      ## Get the associated origin mean for each site
    p_lad_sites_meanOrigins <- p_list$origins[site_origins_id]
    
      ## Compute the hierarchical log-likelihood
    log_likelihood_sites <- dnorm(p_list$sites, mean = p_lad_sites_meanOrigins, sd = p_list$sigma_sites, log = TRUE)
    
    
    
    # 5. Compute the data log-likelihood
    log_likelihood_data <- dnorm(residuals, mean = 0, sd = p_list$sigma, log = TRUE) 
    
    
    
    # 6. Compute the total log likelihood
    
    # If pointwise log-likelihood (used for WAIC), return the sum of the log likelihood ONLY from observations
    # For WAIC, the pointwise log-likelihoods should be based only on the likelihood of observed data, 
    # not including the random effects prior (i.e., their contribution to the joint likelihood). 
    # This is consistent with the fact that WAIC is estimating expected predictive performance — 
    # you want to assess how well the model predicts new observed data, not how well it fits the prior on latent variables.

    if (pointwise) {
      return(list(
        "residuals" = residuals,
        "loglikelihood" = log_likelihood_data
      ))
    }
    
    
    # Else, during the MCMC sampling, return the log-posterior scalar
    # As the sum of the log-likelihood from observed data and from random hierarchical effects
    log_posterior <- sum(log_likelihood_data) + sum(log_likelihood_origins) + sum(log_likelihood_sites)
    
    return(log_posterior) 
  }
  
  
  ## FUNCTION TO CREATE THE PRIORS ----
  # First parameters are the half cauchy priors (used to model SD)
  # Others are uniform priors
  createCombinedPriors <- function(halfcauchy_scale,
                                   uniform_lower, uniform_upper) {
    
    # Some tests
    if (length(uniform_lower) != length(uniform_upper)) stop("Not same length of lower and upper uniform bounds")
    if (any(uniform_upper < uniform_lower)) stop("Some lower bounds are greater than upper uniform bounds")
    
    # Get number of priors for both half cauchy and uniforms
    n_params_halfcauchy <- length(halfcauchy_scale)
    n_params_uniform <- length(uniform_upper)
    
    # Log-prior density function
    prior_density <- function(params) {
      
      # Extract parameters with a ...
      
      ## Half-Cauchy prior (i.e. the first n_params_halfcauchy parameters)
      halfcauchy_params <- params[1:n_params_halfcauchy]
      
      ## Uniform prior (the other parameters)
      uniform_params <- params[(n_params_halfcauchy+1):length(params)]
      
      
      # Check parameters for ...
      
      ## Half-Cauchy scale must be positive
      if (any(halfcauchy_params <= 0)) return(-Inf)
      
      ## Uniform is bounded by lower and upper
      if (any(uniform_params < uniform_lower) || any(uniform_params > uniform_upper)) {
        return(-Inf)
      }
      
      
      # Log-density functions for ...
      
      ## Half-Cauchy priors
      log_prior_halfcauchy <- sum(extraDistr::dhcauchy(
        halfcauchy_params,
        sigma = halfcauchy_scale,
        log = TRUE
      ))
      
      ## Uniform priors
      
      # Uniform priors within bounded support do not need an explicit density term in most Bayesian sampling frameworks like BayesianTools.
      # When using a uniform prior within known bounds, the prior density is constant (i.e., flat).
      # In Bayesian inference, adding a constant to the log-posterior doesn’t affect the MCMC sampling, 
      # because MCMC uses ratios of densities, and constants cancel out.
      # Therefore, as long as we enforce that the uniform parameters stay within bounds, and
      # we return -Inf (log(0)) if they fall outside the bounds
      # we don’t need to explicitly add a log-density term for uniform priors
      
      
      # Total log-prior (which is thus only the log prior of the halfCauchy distributions)
      return(log_prior_halfcauchy)
    }
    
    
    # Prior sampler function
    prior_sampler <- function(n=1) {
      
      # Sample from Half-Cauchy
      halfcauchy_samples <- sapply(
        1:n_params_halfcauchy, 
        function(i) return(extraDistr::rhcauchy(n, sigma = halfcauchy_scale[i]))
      )
      
      # Sample from Uniform
      uniform_samples <- sapply(
        1:n_params_uniform, 
        function(i) return(runif(n, min = uniform_lower[i], max = uniform_upper[i]))
      )
      
      # Combine
      if (n>1) {
        samples <- t(cbind(halfcauchy_samples, uniform_samples))
      } else {
        samples <- c(halfcauchy_samples, uniform_samples)
      }
      return(samples)
    }
    
    # Combine lower and upper bounds of both halfCauchy and uniform priors
    lower <- c(rep(0, n_params_halfcauchy), uniform_lower)
    upper <- c(rep(Inf, n_params_halfcauchy), uniform_upper)
    best <- c(halfcauchy_scale, (uniform_upper + uniform_lower) / 2)
    
    return(BayesianTools::createPrior(
      density = prior_density,
      sampler = prior_sampler,
      lower = lower,
      upper = upper,
      best = best
    ))
  }
  
  
  
  
  # SCRIPT ----
  
  ## Filter sensors ----
  
  # Remove sensors that did not converged
  for (site_name in names(data_sensors)) {
    
    data_sensors[[site_name]] <- data_sensors[[site_name]] %>%
      dplyr::left_join(
        output_lad_method1,
        by = c("plot" = "site", "id" = "id_sensor")
      ) %>%
      dplyr::filter(converged)
    
    data_calib[[site_name]]$sensors <- data_calib[[site_name]]$sensors %>%
      dplyr::left_join(
        output_lad_method1 %>% dplyr::filter(site == site_name),
        by = "id_sensor"
      ) %>%
      dplyr::filter(converged)
    
  }
  
  message(
    paste0(
      sum(!output_lad_method1$converged), "/",
      nrow(output_lad_method1), " sensors had been removed"
    )
  )
  
  
  # Remove sites without sensors
  sites_nosensors <- unlist(data_sensors %>% purrr::map(~nrow(.x) == 0))
  data_sensors <- data_sensors[!sites_nosensors]
  
  message(
    paste0(
      sum(sites_nosensors), " sites had been removed : ",
      paste(names(sites_nosensors)[sites_nosensors], collapse = ", ")
    )
  )
  
  # Filter data plots
  data_plots <- data_plots %>% 
    dplyr::filter(name %in% names(data_sensors))
  
  
  
  
  
  ## Base variables ----
  origin_names <- unique(data_plots$origin)
  site_names <- data_plots$name
  
  site_origins <- data_plots$origin # Origins of all the site
  site_origins_id <- match(site_origins, origin_names) # Corresponding id in origin_names for all the sites
  
  n_species <- length(species2calib)
  n_origins <- length(origin_names)
  n_sites <- length(site_names)
  
  
  
  
  
  ## Define parameters and priors ----
  # HALFCAUCHY PRIORS MUST BE THE FIRST ONES (typically for modelling SD)
  
  # BASE PARAMETER
  par_names <- c(
    "sigma" # model sigma SD
  )
  
  # RANDOM EFFECT
  
  ## SD of the random origin effect (hyperparameter)
  if (mod_design$origin_rd_effect) {
    par_names <- c(par_names, "sigma_origin")
  }
  
  ## SD of the random site effect (hyperparameter)
  if (mod_design$sites_rd_effect) {
    par_names <- c(par_names, "sigma_site")
  }
  
  ## Random origins intercept
  if (mod_design$origin_rd_effect) {
    par_names <- c(par_names, paste0("origin.", origin_names))
  }
  
  ## Random sites intercept
  if (mod_design$sites_rd_effect) {
    par_names <- c(par_names, paste0("site.", site_names))
  }
  
  prior_halfcauchy_S <- rep(5, times = 1 + sum(mod_design$origin_rd_effect, mod_design$sites_rd_effect)) # SD parameters
  
  prior_uniform_LB <- rep(-10, times = n_origins*mod_design$origin_rd_effect + n_sites * mod_design$sites_rd_effect) # site and origin random effects
  prior_uniform_UB <- rep(10, n_origins*mod_design$origin_rd_effect + n_sites * mod_design$sites_rd_effect) # site and origin random effects
  
  
  
  # PARAMS DEPENDING ON MODEL DEFINITION
  
  if (mod_design$species_specific) {
    
    # species-specific intercept
    par_names <- c(par_names, paste0("species.", species2calib)) 
    prior_uniform_LB <- c(prior_uniform_LB, rep(0.1, times=n_species))
    prior_uniform_UB <- c(prior_uniform_UB, rep(5, times=n_species))
    
  } else {
    
    # a single intercept for all species
    par_names <- c(par_names, "intercept")
    prior_uniform_LB <- c(prior_uniform_LB, 0.1)
    prior_uniform_UB <- c(prior_uniform_UB, 5)
    
  }
  
  
  if (mod_design$dbh_effect) {
    
    if (mod_design$species_specific) {
      
      # dbh effect for each species
      par_names <- c(par_names, paste0("dbh.", species2calib)) 
      prior_uniform_LB <- c(prior_uniform_LB, rep(-0.01, times=n_species))
      prior_uniform_UB <- c(prior_uniform_UB, rep(0.01, times=n_species))
      
    } else {
      
      # a single dbh effect for all species
      par_names <- c(par_names, "dbh")
      prior_uniform_LB <- c(prior_uniform_LB, -0.01)
      prior_uniform_UB <- c(prior_uniform_UB, 0.01)
      
    }
    
  } # Otherwise, no dbh effect
  
  
  if (mod_design$compet_effect) {
    if (mod_design$species_specific) {
      # competition effect for each species
      par_names <- c(par_names, paste0("compet.", species2calib)) 
      prior_uniform_LB <- c(prior_uniform_LB, rep(-0.01, times=n_species))
      prior_uniform_UB <- c(prior_uniform_UB, rep(0.01, times=n_species))
    } else {
      # a single competition effect for all species
      par_names <- c(par_names, "compet")
      prior_uniform_LB <- c(prior_uniform_LB, -0.01)
      prior_uniform_UB <- c(prior_uniform_UB, 0.01)
    }
  } # Otherwise, no competition effect
  
  
  
  
  ## Create priors ----
  prior <- createCombinedPriors(halfcauchy_scale = prior_halfcauchy_S,
                                uniform_lower = prior_uniform_LB,
                                uniform_upper = prior_uniform_UB)
  
  
  
  
  
  ## Bayesian setup ----
  bayesianSetup <- BayesianTools::createBayesianSetup(compute_log_posterior, 
                                                      prior, 
                                                      names = par_names)
  
  
  
  ## Return the environment ----
  return(rlang::current_env())
}
