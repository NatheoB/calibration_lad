initialise_models <- function(exp_design,
                              data_sensors,
                              data_plots,
                              data_calib,
                              data_rad,
                              output_lad_method1,
                              species2calib) {
  
  
  ## Create the Bayesian setup for each model ----
  ids_simu <- exp_design$id_simu
  model_setups <- setNames(vector("list", length(ids_simu)), ids_simu)
  for (i in 1:length(ids_simu)) {
    
    id_simu <- ids_simu[i]
    
    model_setups[[id_simu]] <- initialise_model(exp_design[i,],
                                                data_sensors, 
                                                data_plots,
                                                data_calib, 
                                                data_rad,
                                                output_lad_method1,
                                                species2calib)
  }
  
  
  ## Return the model setups ----
  return(model_setups)
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
                               lad_species
                               # lad_dbh,
                               # lad_compet
                               ) {
    
    # Species-specific parameters ----
    # Here the order of the species are really important
    # It is the same as in the species2calib and in p_lad_intercepts
    sp_intercepts <- data.frame(
      species_calib = species2calib,
      intercept = lad_species
      # dbh_coef = lad_dbh,
      # compet_coef = lad_compet
    )
    
    # Random site effect ----
    # Here, the order is really important, and is the same between site_names and lad_sites vectors
    site_intercept <- lad_sites[which(site_names == site)]
    
    # Compute the tree LAD ----
    # Here, only add the site random effect, as the origin effect in already included in the estimation of the site effect (hierarchical structure)
    tmp_trees <- data_calib[[site]]$trees %>% 
      dplyr::left_join(sp_intercepts, by = "species_calib") %>% 
      dplyr::mutate( 
        eta = site_intercept + intercept,
        crown_lad = exp(eta)
      )
    
    # Get the plot info ----
    tmp_plot <- data_plots %>% 
      dplyr::filter(name == site)
    
    # Run SamsaraLight ----
    tmp_out_samsalight <- 
      sl_run(tmp_trees, 
             data_rad[[site]],
             sensors = data_calib[[site]]$sensors, 
             sensors_only = TRUE,
             latitude = tmp_plot$latitude, 
             slope = tmp_plot$slope, 
             aspect = tmp_plot$aspect, 
             north_to_x_cw = tmp_plot$northToX,
             start_day = 1, 
             end_day = 365,
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
    
    # Return sl sensors output for the plot ----
    tmp_out_samsalight$output$sensors
    
  }
  
  ## FUNCTION TO RUN SAMSARALIGHT ON ALL STANDS ----
  compute_pacl_residuals <- function(data_sensors,
                                     data_plots,
                                     data_calib,
                                     data_rad,
                                     lad_sites,
                                     lad_species,
                                     # lad_dbh,
                                     # lad_compet,
                                     print.pb) {
    
    # Output list
    out_residuals_list <- vector("list", length(site_names))
    
    # Initialize progress bar
    i <- 1
    
    if (print.pb) {
      pb <- txtProgressBar(min = 0, max = length(site_names),
                           style = 3, width = 50, char = "=")
    }
    
    # Run for each site ----
    for (site in site_names) {
      
      ## Sl sensors output for a given LAD model in a given plot ----
      tmp_out_sl <- run_sl_standXlad(
        data_plots,
        data_calib, 
        data_rad,
        site, 
        lad_sites,
        lad_species
        # lad_dbh,
        # lad_compet
      )
      
      # Compute PACL residuals the sensors ----
      out_residuals_list[[i]] <- dplyr::left_join(
        
        # Estimated pacl from virtual sensors
        tmp_out_sl %>%
          dplyr::select(id_sensor, pacl_sl = pacl),
        
        # Measured pacl from field sensor
        data_sensors[[site]] %>%
          dplyr::select(id_sensor = id, pacl_field = PACLtotal),
        
        by = "id_sensor"
      ) %>% 
        
        dplyr::mutate(residuals = pacl_sl - pacl_field,
                      site = site) %>% 
        dplyr::select(site, id_sensor, pacl_sl, pacl_field, residuals)
      
      
      # Update the progress bar
      if (print.pb) setTxtProgressBar(pb, i)
      i <- i+1
    }
    if (print.pb) close(pb)
    
    # Return residuals for all sensors ----
    dplyr::bind_rows(out_residuals_list)
  }
  
  
  
  ## FUNCTION TO CREATE THE SUBVECTOR OF PARAMETERS ----
  split_parameters_vector <- function(p) {
    
    i_param <- 1 # Index to follow the start parameter during the vector splitting
    
    # 1. Standard deviations ----
    
    ## 1.1. Model residuals ----
    p_sigma <- p[i_param]
    i_param <- i_param + 1
    
    ## 1.2. Origin random effect----
    if (mod_design$origin_rd_effect) {
      p_sigma_origin <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_sigma_origin <- 0
    }
    
    ## 1.3. Site random effect ----
    if (mod_design$site_rd_effect) {
      p_sigma_site <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_sigma_site <- 0
    }
    
    
    # 2. Latent variables of random effects ----
    
    ## 2.1. Origin random effect----
    if (mod_design$origin_rd_effect) {
      p_z_origins <- p[i_param + (1:n_origins) - 1]
      i_param <- i_param + n_origins
    } else {
      p_z_origins <- rep(0, n_origins)
    }
    
    ## 2.2. Site random effect ----
    if (mod_design$site_rd_effect) {
      p_z_sites <- p[i_param + (1:n_sites) - 1]
      i_param <- i_param + n_sites
    } else {
      p_z_sites <- rep(0, n_sites)
    }
    
    
    
    # 3. LAD intercept ----
    
    if (mod_design$species_specific) {
      # species-specific intercept
      p_lad_species <- p[i_param + (1:n_species) - 1]
      i_param <- i_param + n_species

    } else {
      # a single intercept for all species
      p_lad_species <- rep(p[i_param], n_species)
      i_param <- i_param + 1
    }
    
    
    # List with all different group of parameters ----
    return(list(
      "sigma_log" = p_sigma,
      "sigma_origin_log" = p_sigma_origin,
      "sigma_site_log" = p_sigma_site,
      "z_origins" = p_z_origins,
      "z_sites" = p_z_sites,
      "species" = p_lad_species
    ))
    
    
      
    # # 2. Random effects
    # 
    # ## 2.1. Standard deviations
    # if (mod_design$origin_rd_effect) {
    #   p_sigma_origins <- p[i_param] # SD of the origin effect
    #   i_param <- i_param + 1
    #   
    # } else {
    #   # Otherwise, NULL standard deviation
    #   p_sigma_origins <- NULL
    # }
    # 
    # 
    # if (mod_design$site_rd_effect) {
    #   p_sigma_sites <- p[i_param] # SD of the site effect
    #   i_param <- i_param + 1
    # 
    # } else {
    #   # Otherwise, NULL standard deviation
    #   p_sigma_sites <- NULL
    # }
    # 
    # 
    # ## 2.2. Random intercepts
    # if (mod_design$origin_rd_effect) {
    #   p_lad_origins <- p[i_param + (1:n_origins) - 1] # Origins intercepts
    #   i_param <- i_param + n_origins
    #   
    # } else {
    #   # Otherwise, random origins effect is always 0
    #   p_lad_origins <- rep(0, n_origins)
    # }
    # 
    # 
    # if (mod_design$site_rd_effect) {
    #   p_lad_sites <- p[i_param + (1:n_sites) - 1] # Site intercepts
    #   i_param <- i_param + n_sites
    #   
    # } else {
    #   # Otherwise, random site effect is always 0
    #   p_lad_sites <- rep(0, n_sites)
    # }
    # 
    
    # # 3. Species-specific parameters
    # 
    # ## 3.1. Intercept
    # 
    # if (mod_design$species_specific) {
    #   # species-specific intercept
    #   p_lad_species <- p[i_param + (1:n_species) - 1]
    #   i_param <- i_param + n_species
    #   
    # } else {
    #   # a single intercept for all species
    #   p_lad_species <- rep(p[i_param], n_species)
    #   i_param <- i_param + 1
    # }
    # 
    # 
    # ## 3.2. DBH effect
    # if (mod_design$dbh_effect) {
    #   
    #   if (mod_design$species_specific) {
    #     # dbh effect for each species
    #     p_lad_dbh <- p[i_param + (1:n_species) - 1]
    #     i_param <- i_param + n_species
    #     
    #   } else {
    #     # a single dbh effect for all species
    #     p_lad_dbh <- rep(p[i_param], n_species)
    #     i_param <- i_param + 1
    #   }
    #   
    # } else {
    #   # Otherwise, dbh effect is always 0
    #   p_lad_dbh <- rep(0, n_species)
    # } 
    # 
    # 
    # if (mod_design$compet_effect) {
    #   
    #   if (mod_design$species_specific) {
    #     # compet effect for each species
    #     p_lad_compet <- p[i_param + (1:n_species) - 1]
    #     i_param <- i_param + n_species
    #     
    #   } else {
    #     # a single compet effect for all species
    #     p_lad_compet <- rep(p[i_param], n_species)
    #     i_param <- i_param + 1
    #   }
    #   
    # } else {
    #   # Otherwise, compet effect is always 0
    #   p_lad_compet <- rep(0, n_species)
    # } 
    # 
    # 
    # # List with all different group of parameters 
    # return(list(
    #   "sigma" = p_sigma,
    #   "sigma_origins" = p_sigma_origins,
    #   "sigma_sites" = p_sigma_sites,
    #   "origins" = p_lad_origins,
    #   "sites" = p_lad_sites,
    #   "species" = p_lad_species,
    #   "dbh" = p_lad_dbh,
    #   "compet" = p_lad_compet
    # ))
  }
  
  
  ## FUNCTION TO COMPUTE LOG LIKELIHOOD OF THE DATA ----
  compute_log_likelihood <- function(p, pointwise = FALSE, print.pb = FALSE) {
    
    # 1. Parameters of the LAD model ---- 
    
    ## 1.1. Unpack the parameter vector ----
    p_list <- split_parameters_vector(p)
  
    ## 1.2. Unlog the sigma parameters ----
    p_sigma <- exp(p_list$sigma_log)
    
    p_sigma_origin <- exp(p_list$sigma_origin_log)
    p_sigma_site <- exp(p_list$sigma_site_log)
    
    ## 1.3. Check for invalid parameter values ----
    # if (any(!is.finite(unlist(p_list)))) return(-Inf)

    
    # 2. Compute the non-centered hierarchical site/origin random effect ----
    
    ## 2.1. Compute the origins mean ----
    mean_origins <- p_sigma_origin * p_list$z_origins  # z_origin ~ N(0,1)
    
    ## 2.2. Get the associated origin mean for each site ----
    mean_sites <- mean_origins[site_origins_id]
    
    ## 2.3. Site effect nested in the origin effect ----
    p_sites <- mean_sites + p_sigma_site * p_list$z_sites  # z_site ~ N(0,1)
    
    
    
    # 3. Compute the residuals with the given lad parameters ----
    
    ## 3.1. get residuals by running SmasaraLight ----
    out_sl <- compute_pacl_residuals(data_sensors,
                                     data_plots,
                                     data_calib, 
                                     data_rad, 
                                     p_sites,
                                     p_list$species,
                                     # p_list$dbh,
                                     # p_list$compet,
                                     print.pb)
    residuals <- out_sl$residuals
    
    ## 3.2. Check for invalid residuals
    # if (any(!is.finite(residuals))) return(-Inf)

    
    
    # 4. Compute the log-likelihood
    
    ## 4.1. Total log-likelihood of the data ----
    # During the MCMC sampling,
    # We combine this log-likelihood with the log-priors computed by BayesianTools from priors defined in the BAyesian Setup 
    log_likelihood_data <- dnorm(residuals, 
                                 mean = 0, sd = p_sigma, 
                                 log = TRUE)
    
    ## 4.2. Check for invalid log-likelihood ----
    # if (any(!is.finite(log_likelihood_data))) return(-Inf)
  
    
    # 5. Return summed or pointwise log-likelihood ----
    
    ## 5.1. Pointwise log-likelihood and residuals ----
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
    
    
    ## 5.2. Otherwise, return the log-likelihood scalar ----
    return(sum(log_likelihood_data))
  }
    
  
  
  ## FUNCTION TO CREATE THE PRIORS ----
  # First parameters are the half cauchy priors (used to model SD)
  # Then we have normal priors
  # Last parameters have uniform priors
  createCombinedPriors <- function(halfcauchy_scale = numeric(0),
                                   normal_mean = numeric(0),
                                   normal_sd = numeric(0),
                                   uniform_lower = numeric(0),
                                   uniform_upper = numeric(0)) {
    
    # Some tests
    if (length(uniform_lower) != length(uniform_upper)) 
      stop("uniform_lower and uniform_upper must have same length")
    if (any(uniform_upper < uniform_lower)) 
      stop("each element of uniform_upper must be greater than uniform_lower")
    if (length(normal_mean) != length(normal_sd))
      stop("normal_mean and normal_sd must have same length")
    
    
    # Get number of priors for half cauchy, normal and uniform
    n_params_halfcauchy <- length(halfcauchy_scale)
    n_params_normal <- length(normal_mean)
    n_params_uniform <- length(uniform_lower)
    

    # Log-prior density function
    prior_density <- function(params) {
      
      # Extract parameters with a ...
      
      ## Half-Cauchy prior (i.e. the first n_params_halfcauchy parameters)
      halfcauchy_params <- numeric(0)
      if (n_params_halfcauchy > 0) halfcauchy_params <- params[1:n_params_halfcauchy] 
      
      ## Normal prior (i.e. n_params_normal parameters after half-cauchy params)
      normal_params <- numeric(0)
      if (n_params_normal > 0) normal_params <- params[1:n_params_normal + n_params_halfcauchy]
      
      ## Uniform prior (i.e. the last parameters)
      uniform_params <- numeric(0)
      if (n_params_uniform > 0) uniform_params <- params[1:n_params_uniform + n_params_halfcauchy + n_params_normal]
      
      
      # Check parameters for ...
      
      ## Half-Cauchy scale must be positive
      if (any(halfcauchy_params <= 0)) return(-Inf)
      
      ## Uniform is bounded by lower and upper
      if (any(uniform_params < uniform_lower) || any(uniform_params > uniform_upper)) return(-Inf)
      
      # Log-density functions for ...
      
      ## Half-Cauchy priors
      log_prior_halfcauchy <- 0
      if (n_params_halfcauchy > 0) {
        log_prior_halfcauchy <- sum(extraDistr::dhcauchy(
          halfcauchy_params,
          sigma = halfcauchy_scale,
          log = TRUE
        ))
      }
        
      
      ## Normal priors
      log_prior_normal <- 0
      if (n_params_normal > 0) {
        log_prior_normal <- sum(dnorm(
          normal_params, 
          mean = normal_mean, 
          sd = normal_sd, 
          log = TRUE
        ))
      }
      
      
      ## Uniform priors
      
      # Uniform priors within bounded support do not need an explicit density term in most Bayesian sampling frameworks like BayesianTools.
      # When using a uniform prior within known bounds, the prior density is constant (i.e., flat).
      # In Bayesian inference, adding a constant to the log-posterior doesn’t affect the MCMC sampling, 
      # because MCMC uses ratios of densities, and constants cancel out.
      # Therefore, as long as we enforce that the uniform parameters stay within bounds, and
      # we return -Inf (log(0)) if they fall outside the bounds
      # we don’t need to explicitly add a log-density term for uniform priors
      
      
      # Total log-prior
      log_prior_total <- log_prior_halfcauchy + log_prior_normal
      
      return(log_prior_total)
    }
    
    
    # Prior sampler function
    prior_sampler <- function(n=1) {
      
      # Sample from Half-Cauchy
      halfcauchy_samples <- NULL
      if (n_params_halfcauchy > 0) {
        halfcauchy_samples <- sapply(
          1:n_params_halfcauchy, 
          function(i) return(extraDistr::rhcauchy(n, sigma = halfcauchy_scale[i]))
        )
      }
        
      
      # Sample from Normal
      normal_samples <- NULL
      if (n_params_normal > 0) {
        normal_samples <- sapply(
          1:n_params_normal, 
          function(i) return(rnorm(n, mean = normal_mean[i], sd = normal_sd[i]))
        )
      }
      
      # Sample from Uniform
      uniform_samples <- NULL
      if (n_params_uniform > 0) {
        uniform_samples <- sapply(
          1:n_params_uniform, 
          function(i) return(runif(n, min = uniform_lower[i], max = uniform_upper[i]))
        )
      }
      

      # Combine samples
      if (n == 1) {
        samples <- c(halfcauchy_samples, normal_samples, uniform_samples)
      } else {
        samples <- cbind(halfcauchy_samples, normal_samples, uniform_samples)
      }
      
      return(samples)
    }
    
    # Create BayesianTools prior object ----
    lower <- c(rep(0, n_params_halfcauchy),
               rep(-Inf, n_params_normal),
               uniform_lower)
    
    upper <- c(rep(Inf, n_params_halfcauchy),
               rep(Inf, n_params_normal),
               uniform_upper)
    
    best <- c(halfcauchy_scale,
              normal_mean,
              (uniform_upper + uniform_lower) / 2)
    
    
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
  
  if (mod_design$filter_sensors) {
    
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
    
  }
  
  
  
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
  # Then normal priors and finally uniform priors
  
  par_halfCauchy <- character(0)
  prior_halfcauchy_S <- numeric(0)
  
  par_normal <- character(0)
  prior_normal_mean <- numeric(0)
  prior_normal_sd <- numeric(0)
  
  par_uniform <- character(0)
  prior_uniform_LB <- numeric(0)
  prior_uniform_UB <- numeric(0)
  
  

  ### Standard deviations ----
  # Here model log of SD with a Normal prior, to ensure for positivity
  # Faster and easier to converge than half cauchy prior
  # Especially with nested random effects
  
  #### Model residuals ----
  par_normal <- c(par_normal, "sigma_log")
  prior_normal_mean <- c(prior_normal_mean, log(0.3))
  prior_normal_sd <- c(prior_normal_sd, 0.5)
  
  #### Origin random effect ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, "sigma_origin_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.5)
  }
  
  #### Site random effect ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, "sigma_site_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.5)
  }
  
  # #### Model residuals
  # par_halfCauchy <- c(par_halfCauchy, "sigma")
  # prior_halfcauchy_S <- c(prior_halfcauchy_S, 3)
  # 
  # #### Origin random effect
  # if (mod_design$origin_rd_effect) {
  #   par_halfCauchy <- c(par_halfCauchy, "sigma_origin")
  #   prior_halfcauchy_S <- c(prior_halfcauchy_S, 5)
  # }
  # 
  # #### Site random effect
  # if (mod_design$site_rd_effect) {
  #   par_halfCauchy <- c(par_halfCauchy, "sigma_site")
  #   prior_halfcauchy_S <- c(prior_halfcauchy_S, 5)
  # }
 
  
  ### Latent variables of random effects ----
  # Latent variables N(0,1)
  
  #### Origin random effect ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, paste0("origin.", origin_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_origins))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_origins))
  }
  
  #### Site random effect ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, paste0("site.", site_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_sites))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_sites))
  }
  
  
  ### LAD intercept ----
  # (based on the best LAD guess from the methodology 1)
  # Here, we estimate the log(LAD), to ensure for LAD positivity
  lad_guess_method1 <- mean(output_lad_method1$best_lad[output_lad_method1$converged])
  
  if (lad_guess_method1 <= 0) stop("best lad guess from methodo is lower or equal to 0")
  lad_guess_method1_log <- log(lad_guess_method1)
  
  if (mod_design$species_specific) {
    
    # species-specific intercept
    par_normal <- c(par_normal, paste0("species.", species2calib)) 
    prior_normal_mean <- c(prior_normal_mean, rep(lad_guess_method1_log, times=n_species))
    prior_normal_sd <- c(prior_normal_sd, rep(0.5, times=n_species))
    
  } else {
    
    # a single intercept for all species
    par_normal <- c(par_normal, "intercept")
    prior_normal_mean <- c(prior_normal_mean, lad_guess_method1_log)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
    
  }
  
  
  
  # DBH effect
  if (mod_design$dbh_effect) {
    stop("-- DEBUG -- : DBH effect does not work")
    
    # if (mod_design$species_specific) {
    #   
    #   # dbh effect for each species
    #   par_names <- c(par_names, paste0("dbh.", species2calib)) 
    #   prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
    #   prior_normal_sd <- c(prior_normal_sd, rep(0.1, times=n_species))
    #   
    # } else {
    #   
    #   # a single dbh effect for all species
    #   par_names <- c(par_names, "dbh")
    #   prior_normal_mean <- c(prior_normal_mean, 0)
    #   prior_normal_sd <- c(prior_normal_sd, 0.1)
    #   
    # }
    
  } # Otherwise, no dbh effect

  
  ## Competition effect
  if (mod_design$compet_effect) {
    stop("DEBUG: competition effect do not work now")
    
    # if (mod_design$species_specific) {
    #   # competition effect for each species
    #   par_names <- c(par_names, paste0("compet.", species2calib)) 
    #   prior_uniform_LB <- c(prior_uniform_LB, rep(-0.01, times=n_species))
    #   prior_uniform_UB <- c(prior_uniform_UB, rep(0.01, times=n_species))
    # } else {
    #   # a single competition effect for all species
    #   par_names <- c(par_names, "compet")
    #   prior_uniform_LB <- c(prior_uniform_LB, -0.01)
    #   prior_uniform_UB <- c(prior_uniform_UB, 0.01)
    # }
  } # Otherwise, no competition effect
  
  
  ## Combine parameter names ----
  par_names <- c(par_halfCauchy, par_normal, par_uniform)
  
  
  ## Create priors ----
  priors <- createCombinedPriors(
    halfcauchy_scale = prior_halfcauchy_S,
    normal_mean = prior_normal_mean,
    normal_sd = prior_normal_sd,
    uniform_lower = prior_uniform_LB,
    uniform_upper = prior_uniform_UB
  )
  
  
  ## Bayesian setup ----
  bayesianSetup <- BayesianTools::createBayesianSetup(
    likelihood = compute_log_likelihood, 
    prior = priors,
    names = par_names
  )
  
  
  ## Return the environment ----
  return(rlang::current_env())
}
