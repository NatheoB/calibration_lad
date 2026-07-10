initialise_models <- function(exp_design,
                              data_plots,
                              data_stands,
                              data_rad,
                              sensors_punobs,
                              prior_lad,
                              n_threads) {
  
  
  ## Create the Bayesian setup for each model ----
  ids_model <- exp_design$id_model
  model_setups <- setNames(vector("list", length(ids_model)), ids_model)
  for (i in 1:length(ids_model)) {
    
    id_model <- ids_model[i]
    
    model_setups[[id_model]] <- initialise_model(exp_design[i,],
                                                 data_plots,
                                                 data_stands,
                                                 data_rad,
                                                 sensors_punobs,
                                                 prior_lad,
                                                 n_threads)
  }
  
  
  ## Return the model setups ----
  return(model_setups)
}



initialise_model <- function(mod_design,
                             data_plots,
                             data_stands,
                             data_rad,
                             sensors_punobs,
                             prior_lad,
                             n_threads) {
  
  # BayesianTools consider global variables during the MCMC sampling
  # Thus, we need to create an environment with all the functions and variables needed for MCMC sampling
  
  # SUB-FUNCTIONS ----
  
  ## FUNCTION TO RUN SAMSARALIGHT ----
  run_sl_standXlad <- function(data_stands, 
                               data_rad,
                               site, 
                               lad_random_per_sites,
                               lad_intercept_per_sp,
                               lad_shadetol_per_sp,
                               lad_dbh_per_sp,
                               lad_compet_per_sp,
                               lad_shadetolXdbh_per_sp,
                               lad_shadetolXcompet_per_sp,
                               lad_dbhXcompet_per_sp,
                               lad_shadetolXdbhXcompet_per_sp
                               ) {
    
    # Species-specific parameters ----
    # Here the order of the species modalities is really important
    # It is the same as in the species2calib, in lad_intercept_per_species...
    # And also the same in is_species_gymno when we added gymnosperm coefficients to species-specific parameters
    sp_coefs <- data.frame(
      species = species2calib,
      intercept = lad_intercept_per_sp,
      beta_shadetol = lad_shadetol_per_sp,
      beta_dbh = lad_dbh_per_sp,
      beta_compet = lad_compet_per_sp,
      beta_shadetolXdbh = lad_shadetolXdbh_per_sp,
      beta_shadetolXcompet = lad_shadetolXcompet_per_sp,
      beta_dbhXcompet = lad_dbhXcompet_per_sp,
      beta_shadetolXdbhXcompet = lad_shadetolXdbhXcompet_per_sp
    )

    
    # Random site effect ----
    # Here, the order is really important, and is the same between site_names and lad_sites vectors
    random_site <- lad_random_per_sites[which(site_names == site)]
    
    # Compute the tree LAD ----
    # Here, only add the site random effect, as the origin effect in already included in the estimation of the site effect (hierarchical structure)
    # Same, we did also consider the gymnosperm effect, adding them to each species-specific coefficient
    tmp_stand <- data_stands[[site]]
    
    tmp_stand$trees <- tmp_stand$trees %>% 
      dplyr::left_join(sp_coefs, by = "species") %>% 
      dplyr::mutate( 
        eta = random_site + 
          intercept + 
          beta_shadetol * shadetolstd + 
          beta_dbh * dbhstd + 
          beta_compet * competstd + 
          beta_shadetolXdbh * shadetolstd * dbhstd + 
          beta_shadetolXcompet * shadetolstd * competstd + 
          beta_dbhXcompet * dbhstd * competstd + 
          beta_shadetolXdbhXcompet * shadetolstd * dbhstd * competstd,
        
        crown_lad = exp(eta)
      )
    
    
    # Run SamsaraLight ----
    tmp_out_samsalight <- SamsaRaLight::run_sl(
      tmp_stand, 
      data_rad[[site]],
      sensors_only = TRUE,
      detailed_output = FALSE,
      parallel_mode = TRUE,
      n_threads = n_threads,
      verbose = FALSE
      )
    
    # Return sl sensors output for the plot ----
    tmp_out_samsalight$output$light$sensors
    
  }
  
  ## FUNCTION TO RUN SAMSARALIGHT ON ALL STANDS ----
  compute_pacl_residuals <- function(data_stands,
                                     data_rad,
                                     lad_random_per_sites,
                                     lad_intercept_per_sp,
                                     lad_shadetol_per_sp,
                                     lad_dbh_per_sp,
                                     lad_compet_per_sp,
                                     lad_shadetolXdbh_per_sp,
                                     lad_shadetolXcompet_per_sp,
                                     lad_dbhXcompet_per_sp,
                                     lad_shadetolXdbhXcompet_per_sp,
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
        data_stands,
        data_rad,
        site, 
        lad_random_per_sites,
        lad_intercept_per_sp,
        lad_shadetol_per_sp,
        lad_dbh_per_sp,
        lad_compet_per_sp,
        lad_shadetolXdbh_per_sp,
        lad_shadetolXcompet_per_sp,
        lad_dbhXcompet_per_sp,
        lad_shadetolXdbhXcompet_per_sp
      )
      
      # Compute PACL residuals the sensors ----
      out_residuals_list[[i]] <- dplyr::left_join(
        
        # Estimated pacl from virtual sensors
        tmp_out_sl %>%
          dplyr::select(id_sensor, pacl_sl = pacl),
        
        # Measured pacl from field sensor
        data_stands[[site]]$sensors %>%
          dplyr::select(id_sensor, pacl_field = PACLtotal),
        
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
    p_sigma_mod <- p[i_param]
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
    
    ## 1.4. Species random effect ----
    if (mod_design$species_rd_effect) {
      p_sigma_sp <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_sigma_sp <- 0
    }
    
    
    # 2. Latent variables for random effects ----
    
    ## 2.1. Origins ----
    if (mod_design$origin_rd_effect) {
      p_z_origins <- p[i_param + (1:n_origins) - 1]
      i_param <- i_param + n_origins
    } else {
      p_z_origins <- rep(0, n_origins)
    }
    
    ## 2.2. Sites ----
    if (mod_design$site_rd_effect) {
      p_z_sites <- p[i_param + (1:n_sites) - 1]
      i_param <- i_param + n_sites
    } else {
      p_z_sites <- rep(0, n_sites)
    }
    
    ## 2.3. Species ----
    if (mod_design$species_rd_effect) {
      p_z_species <- p[i_param + (1:n_species) - 1]
      i_param <- i_param + n_species
    } else {
      p_z_species <- rep(0, n_species)
    }
    
    
    # 3. Predictors effect ----
    
    ## 3.1. Intercept ----
    if (!is.na(mod_design$grouping_var)) {
      
      p_intercept_groups <- p[i_param + (1:n_groups) - 1]
      i_param <- i_param + n_groups
      
    } else {
      
      p_intercept_groups <- p[i_param]
      i_param <- i_param + 1
    }
    
    ## 3.2. Slopes ----
    
    ### 3.2.1. SHADETOL effect ----
    if (mod_design$shadetol_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_shadetol_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_shadetol_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no dbh effect
      p_shadetol_groups <- rep(0, times = n_groups)
    }
    
    
    ### 3.2.2. DBH effect ----
    if (mod_design$dbh_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_dbh_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_dbh_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no dbh effect
      p_dbh_groups <- rep(0, times = n_groups)
    }
    
    
    
    ### 3.2.3 Competition effect ----
    if (mod_design$compet_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_compet_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_compet_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no compet effect
      p_compet_groups <- rep(0, times = n_groups)
    }
    
    
    ## 3.3. Interactions ----
    
    ### 3.3.1. SHADETOL X DBH interaction ----
    if (mod_design$shadetol_effect & mod_design$dbh_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_shadetolXdbh_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_shadetolXdbh_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no shadetol X dbh effect
      p_shadetolXdbh_groups <- rep(0, times = n_groups)
    }
    
    
    ### 3.3.2. SHADETOL X COMPET interaction ----
    if (mod_design$shadetol_effect & mod_design$compet_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_shadetolXcompet_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_shadetolXcompet_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no shadetol X compet effect
      p_shadetolXcompet_groups <- rep(0, times = n_groups)
    }
    
    
    ### 3.3.3. DBH X COMPET interaction ----
    if (mod_design$dbh_effect & mod_design$compet_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_dbhXcompet_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        

      } else {
        
        p_dbhXcompet_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no dbh X compet effect
      p_dbhXcompet_groups <- rep(0, times = n_groups)
    }
    
    
    ### 3.3.4. SHADETOL X DBH X COMPET interaction ----
    if (mod_design$shadetol_effect & mod_design$dbh_effect & mod_design$compet_effect) {
      
      if (!is.na(mod_design$grouping_var)) {
        
        p_shadetolXdbhXcompet_groups <- p[i_param + (1:n_groups) - 1]
        i_param <- i_param + n_groups
        
      } else {
        
        p_shadetolXdbhXcompet_groups <- p[i_param]
        i_param <- i_param + 1
      }
      
      
    } else {
      
      # Otherwise, no dbh X compet effect
      p_shadetolXdbhXcompet_groups <- rep(0, times = n_groups)
    }
    
    
    # List with all different group of parameters ----
    return(list(
      
      # Model SD residuals
      "sigma_mod_log" = p_sigma_mod,
      
      # Hierarchical random effect
      "sigma_origin_log" = p_sigma_origin,
      "sigma_site_log" = p_sigma_site,
      "z_random_per_origin" = p_z_origins,
      "z_random_per_site" = p_z_sites,
      
      # Species random effect
      "sigma_sp_log" = p_sigma_sp,
      "z_random_per_sp" = p_z_species,
      
      # Intercept and slope effects
      "intercept_per_group" = p_intercept_groups,
      "shadetol_per_group" = p_shadetol_groups,
      "dbh_per_group" = p_dbh_groups,
      "compet_per_group" = p_compet_groups,
      
      # Interactions
      "shadetolXdbh_per_group" = p_shadetolXdbh_groups,
      "shadetolXcompet_per_group" = p_shadetolXcompet_groups,
      "dbhXcompet_per_group" = p_dbhXcompet_groups,
      "shadetolXdbhXcompet_per_group" = p_shadetolXdbhXcompet_groups
      
    ))
  }
  
  
  ## FUNCTION TO COMPUTE LOG LIKELIHOOD OF THE DATA ----
  compute_log_likelihood <- function(p, 
                                     pointwise = FALSE, 
                                     print.pb = FALSE) {
    
    # 1. Parameters of the LAD model ---- 
    
    ## 1.1. Unpack the parameter vector ----
    p_list <- split_parameters_vector(p)
    
    ## 1.2. Check for invalid parameter values ----
    # if (any(!is.finite(unlist(p_list)))) return(-Inf)
    
    
    
    # 2. Compute the hierarchical site/origin random effects ----
    # non-centered parameterization
    # And within-site effetc of total basal area
    
    ### 2.1. Unlog sigma parameters ----
    p_sigma_origin <- exp(p_list$sigma_origin_log)
    p_sigma_site <- exp(p_list$sigma_site_log)
    
    ## 2.2. Compute the origins mean ----
    mean_origins <- p_sigma_origin * p_list$z_random_per_origin  # z_origin ~ N(0,1)
    
    ## 2.3. Get the associated origin mean for each site ----
    mean_origins_per_site <- mean_origins[id_origin_per_site]
    
    ## 2.4. Site effect nested in the origin effect ----
    p_random_per_site <- mean_origins_per_site + p_sigma_site * p_list$z_random_per_site  # z_site ~ N(0,1)
    
    
    # 3. Compute species random effect ----
    
    ## 3.1. Unlog sigma parameters ----
    p_sigma_sp <- exp(p_list$sigma_sp_log)
    
    ## 3.2. Compute species random effect ----
    p_random_per_sp <- p_sigma_sp * p_list$z_random_per_sp  # z_sp ~ N(0,1)
    
    
    # 4. Get the group-specific intercept and slope effects for each species ----
    p_intercept_per_sp <- p_list$intercept_per_group[id_group_per_species]
    
    p_shadetol_per_sp <- p_list$shadetol_per_group[id_group_per_species]
    p_dbh_per_sp <- p_list$dbh_per_group[id_group_per_species]
    p_compet_per_sp <- p_list$compet_per_group[id_group_per_species]
    
    p_shadetolXdbh_per_sp <- p_list$shadetolXdbh_per_group[id_group_per_species]
    p_shadetolXcompet_per_sp <- p_list$shadetolXcompet_per_group[id_group_per_species]
    p_dbhXcompet_per_sp <- p_list$dbhXcompet_per_group[id_group_per_species]
    p_shadetolXdbhXcompet_per_sp <- p_list$shadetolXdbhXcompet_per_group[id_group_per_species]

    
    # 5. Compute the residuals with the given lad parameters ----
    # (difference between estimated/observed PACL for each sensor)
    
    ## 5.1. get residuals by running SamsaraLight ----
    out_sl <- compute_pacl_residuals(data_stands, 
                                     data_rad, 
                                     p_random_per_site,
                                     p_intercept_per_sp,
                                     p_shadetol_per_sp,
                                     p_dbh_per_sp,
                                     p_compet_per_sp,
                                     p_shadetolXdbh_per_sp,
                                     p_shadetolXcompet_per_sp,
                                     p_dbhXcompet_per_sp,
                                     p_shadetolXdbhXcompet_per_sp,
                                     print.pb)
    residuals <- out_sl$residuals
    
    ## 5.2. Check for invalid residuals ----
    # if (any(!is.finite(residuals))) return(-Inf)
    
    
    
    # 6. Compute the log-likelihood ----
    
    ## 6.1. Unlog sigma parameter ----
    p_sigma_mod <- exp(p_list$sigma_mod_log)
    
    
    ## 6.2. Total log-likelihood of the data ----
    # During the MCMC sampling,
    # We combine this log-likelihood with the log-priors computed by BayesianTools from priors defined in the BAyesian Setup 
    log_likelihood_data <- dnorm(residuals, 
                                 mean = 0, sd = p_sigma_mod, 
                                 log = TRUE)
    
    ## 6.3. Check for invalid log-likelihood ----
    # if (any(!is.finite(log_likelihood_data))) return(-Inf)
    
    
    # 7. Return the output residuals and log-likelihood matrices ----
    
    ## 7.1. Pointwise log-likelihood and residuals ----
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
    
    ## 7.2. Otherwise, return the vector of log-likelihood (i.e. for each sensor) ----
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

  if (mod_design$filter_sensors) {
    
    for (site_name in names(data_stands)) {
      
      data_stands[[site_name]]$sensors <- data_stands[[site_name]]$sensors %>%
        dplyr::left_join(
          sensors_punobs %>% 
            dplyr::filter(site == site_name) %>%
            dplyr::ungroup() %>% 
            dplyr::select(id_sensor, keep),
          by = "id_sensor"
        ) %>%
        dplyr::filter(keep)
      
    }
    
    message(
      paste0(
        sum(!sensors_punobs$keep), "/",
        nrow(sensors_punobs), " sensors had been removed"
      )
    )
    
    
    # Remove sites without sensors
    sites_nosensors <- unlist(data_stands %>% purrr::map(~nrow(.x$sensors) == 0))
    
    data_plots <- data_plots %>% dplyr::filter(!sites_nosensors)
    data_stands <- data_stands[!sites_nosensors]
    
    message(
      paste0(
        sum(sites_nosensors), " sites had been removed : ",
        paste(names(sites_nosensors)[sites_nosensors], collapse = ", ")
      )
    )
    
  }
  
  
  
  ## Standardized variables ----
  
  ### SHADETOL ----
  shadetol_vect <- data_stands %>% 
    purrr::map(~.x$trees) %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(species, shadetol) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(shadetol)
  
  shadetol_mean <- mean(shadetol_vect)
  shadetol_sd <- sd(shadetol_vect)
  
  data_stands <- data_stands %>% 
    purrr::map(~{
      .x$trees <- .x$trees %>% 
        dplyr::mutate(shadetolstd = (shadetol - shadetol_mean) / shadetol_sd)
      .x
    })
  
  ### DBH ----
  dbh_vect <- data_stands %>% 
    purrr::map(~.x$trees$dbh_cm) %>% 
    purrr::reduce(c)
  
  dbh_mean <- mean(dbh_vect)
  dbh_sd <- sd(dbh_vect)
  
  data_stands <- data_stands %>% 
    purrr::map(~{
      .x$trees <- .x$trees %>% 
        dplyr::mutate(dbhstd = (dbh_cm - dbh_mean) / dbh_sd)
      .x
    })
  
  
  ### COMPET ----
  
  # If specified compet effect and variable
  if (mod_design$compet_effect & !is.na(mod_design$compet_var)) {
  
    compet_vect <- data_stands %>% 
      purrr::map(~.x$trees[,mod_design$compet_var]) %>% 
      purrr::reduce(c)
    
    compet_mean <- mean(compet_vect)
    compet_sd <- sd(compet_vect)
    
    data_stands <- data_stands %>% 
      purrr::map(~{
        .x$trees <- .x$trees %>% 
          dplyr::mutate(competstd := ( (!!sym(mod_design$compet_var) - compet_mean) / compet_sd) )
        .x
      })
  
  } else {
    
    # Otherwise, set it to 0 
    data_stands <- data_stands %>% 
      purrr::map(~{
        .x$trees <- .x$trees %>% 
          dplyr::mutate(competstd = 0)
        .x
      })
    
  }
  
  
  ## Global variables ----
  
  ### Site/Origin effect ----
  origin_names <- unique(data_plots$origin)
  site_names <- data_plots$name
  
  site_origins <- data_plots$origin # Origins of all the site
  id_origin_per_site <- match(site_origins, origin_names) # Corresponding id in origin_names for all the sites
  
  n_origins <- length(origin_names)
  n_sites <- length(site_names)
  
  
  ### Species effect ----
  species2calib <- data_stands %>%
    purrr::map(~unique(.x$trees$species)) %>%
    unlist() %>%
    unname() %>%
    unique()
  
  n_species <- length(species2calib)
  
  
  ### Group effect ----
  group_per_species <- NULL
  groups2calib <- NULL
  
  id_group_per_species <- rep(1, times = n_species)
  n_groups <- 1
  
  if (!is.na(mod_design$grouping_var)) {
    
    group_per_species <- data_stands %>%
      
      # Get the grouping variable modality per species
      purrr::map(~.x$trees %>% 
                   dplyr::select(species, group = !!sym(mod_design$grouping_var)) %>% 
                   dplyr::distinct()) %>% 
      dplyr::bind_rows() %>% 
      dplyr::distinct() %>% 
      
      # Order species with the same order as the species vector
      dplyr::mutate(species = factor(species, levels = species2calib)) %>%
      dplyr::arrange() %>%
      dplyr::pull(group)

    groups2calib <- unique(group_per_species)
    id_group_per_species <- match(group_per_species, groups2calib)
    n_groups <- length(groups2calib)
  }

  
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
  par_normal <- c(par_normal, "sigma_mod_log")
  prior_normal_mean <- c(prior_normal_mean, log(0.3))
  prior_normal_sd <- c(prior_normal_sd, 0.3)
  
  
  #### Origin random effect ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, "sigma_origin_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.3)
  }
  
  
  #### Site random effect ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, "sigma_site_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.3)
  }
  
  
  #### Species random effect ----
  if (mod_design$species_rd_effect) {
    par_normal <- c(par_normal, "sigma_sp_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.3)
  }
  
  
  ### Latent variables for random effect ----
  # Non-centered parameterization with latent variables N(0,1)
  # To increase computation time and favour mixing
  
  #### Origins ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, paste0("z_origin.", origin_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_origins))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_origins))
  }
  
  #### Sites ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, paste0("z_site.", site_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_sites))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_sites))
  }

  
  #### Species ----
  if (mod_design$species_rd_effect) {
    par_normal <- c(par_normal, paste0("z_sp.", species2calib))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_species))
  }
  
  
  ### Intercept effect ----
  if (prior_lad <= 0) stop("prior lad is lower or equal to 0")
  prior_lad_log <- log(prior_lad)   # Here, we estimate the log(LAD), to ensure for LAD positivity
  
  if (!is.na(mod_design$grouping_var)) {
    
    # For each group
    par_normal <- c(par_normal, paste0("intercept.", groups2calib))
    prior_normal_mean <- c(prior_normal_mean, rep(prior_lad_log, times = length(groups2calib)))
    prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
    
  } else {
    
    # Mean species
    par_normal <- c(par_normal, "intercept")
    prior_normal_mean <- c(prior_normal_mean, prior_lad_log)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
  }
  
  
  ### SHADETOL effect ----
  # BE CAREFUL : SHADETOL in this model is standardized
  if (mod_design$shadetol_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("shadetol.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "shadetol")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  } # Otherwise, no shadetol effect
  
  
  ### DBH effect ----
  # BE CAREFUL : DBH in this model is standardized
  if (mod_design$dbh_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("dbh.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "dbh")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  } # Otherwise, no dbh effect
  
  
  
  ### COMPET effect ----
  # BE CAREFUL : COMPET variable in this model is standardized
  if (mod_design$compet_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("compet.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "compet")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  } # Otherwise, no compet effect
  
  
  ### Interaction SHADETOL x DBH ----
  if (mod_design$shadetol_effect & mod_design$dbh_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("shadetolXdbh.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "shadetolXdbh")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  }
  
  ### Interaction SHADETOL x COMPET ----
  if (mod_design$shadetol_effect & mod_design$compet_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("shadetolXcompet.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "shadetolXcompet")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  }
  
  ### Interaction DBH X COMPET ----
  if (mod_design$dbh_effect & mod_design$compet_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("dbhXcompet.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "dbhXcompet")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  }
  
  ### Interaction SHADETOL X DBH X COMPET ----
  if (mod_design$shadetol_effect & mod_design$dbh_effect & mod_design$compet_effect) {
    
    if (!is.na(mod_design$grouping_var)) {
      
      # For each group
      par_normal <- c(par_normal, paste0("shadetolXdbhXcompet.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times = length(groups2calib)))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times = length(groups2calib)))
      
    } else {
      
      # Mean species
      par_normal <- c(par_normal, "shadetolXdbhXcompet")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
    }
    
  }
  
  
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
  
  # Initial prior best values
  # Used only for debugging
  best_priors <- setNames(priors$best, par_names)
  
  
  ## Bayesian setup ----
  bayesianSetup <- BayesianTools::createBayesianSetup(
    likelihood = compute_log_likelihood, 
    prior = priors,
    names = par_names
  )
  
  ## Return the environment ----
  return(rlang::current_env())
}
