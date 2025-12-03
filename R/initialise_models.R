initialise_models <- function(exp_design,
                              data_sensors,
                              data_plots,
                              data_calib,
                              data_rad,
                              output_lad_method1,
                              species2calib,
                              parallel = FALSE) {
  
  
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
                             species2calib,
                             parallel = FALSE) {
  
  # BayesianTools consider global variables during the MCMC sampling
  # Thus, we need to create an environment with all the functions and variables needed for MCMC sampling
  
  
  
  # SUB-FUNCTIONS ----
  
  ## FUNCTION TO RUN SAMSARALIGHT ----
  run_sl_standXlad <- function(data_plots,
                               data_calib, 
                               data_rad,
                               site, 
                               lad_random_per_sites,
                               lad_intercept_per_sp,
                               lad_dbh_per_sp
                               # lad_compet
                               ) {
    
    # Species-specific parameters ----
    # Here the order of the species are really important
    # It is the same as in the species2calib and in p_lad_intercepts
    sp_intercepts <- data.frame(
      species_calib = species2calib,
      intercept_sp = lad_intercept_per_sp,
      beta_dbh_sp = lad_dbh_per_sp
      # compet_coef = lad_compet
    )
    
    # Random site effect ----
    # Here, the order is really important, and is the same between site_names and lad_sites vectors
    random_site <- lad_random_per_sites[which(site_names == site)]
    
    # Compute the tree LAD ----
    # Here, only add the site random effect, as the origin effect in already included in the estimation of the site effect (hierarchical structure)
    tmp_trees <- data_calib[[site]]$trees %>% 
      dplyr::left_join(sp_intercepts, by = "species_calib") %>% 
      dplyr::mutate( 
        eta = random_site + intercept_sp + beta_dbh_sp * dbhstd,
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
                                     lad_random_per_sites,
                                     lad_intercept_per_sp,
                                     lad_dbh_per_sp,
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
        lad_random_per_sites,
        lad_intercept_per_sp,
        lad_dbh_per_sp
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
    
    ## 1.4. Intercept effect ----
    ## Hyperparameter SD for intercept 
    if (mod_design$intercept_sp_pooling) {
      p_sigma_intercept <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_sigma_intercept <- 0
    }
    
    ## 1.5. Dbh effect ----
    ## Hyperparameter SD for dbh effect
    if (mod_design$dbh_sp_pooling) {
      p_sigma_dbh <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_sigma_dbh <- 0
    }
    
    
    # 2. Random effects latent variables ----
    
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
    
    
    # 3. Predictors effect ----
    
    ## 3.1. Intercept ----
    
    if (mod_design$intercept_per_sp) {
    
      if (mod_design$intercept_sp_pooling) {
        
        # Hyperparameter mean for species
        ## (i.e. as DBH standardised, mean species LAD at mean DBH)
        p_mean_intercept <- p[i_param]
        i_param <- i_param + 1 
        
        # Species-specific effect for intercept  
        # Latent variables z (for MVN) 
        p_intercept_per_sp <- p[i_param + (1:n_species) - 1]
        i_param <- i_param + n_species
        
      } else {
        
        # otherwise, without partial species spooling
        # fit intercept alpha independently between species
        # Thus set mean to 0 and SD to 1 (i.e. latent variable z become directly the fitted parameter)
        p_mean_intercept <- 0
        p_intercept_per_sp <- p[i_param + (1:n_species) - 1]
        i_param <- i_param + n_species
      }
      
    } else {
      
      # Single intercept over all species
      p_mean_intercept <- 0
      p_intercept_per_sp <- rep(p[i_param], n_species)
      i_param <- i_param + 1
    }
    
    
    
    ## 3.2. DBH effect ----
    if (mod_design$dbh_effect) {
      
      if (mod_design$dbh_per_sp) {
        
        if (mod_design$dbh_sp_pooling) {
          
          # Hyperparameter mean for dbh effect
          p_mean_dbh <- p[i_param]
          i_param <- i_param + 1 
          
          # Species-specific dbh effect
          # Latent variables z (for MVN) 
          p_dbh_per_sp <- p[i_param + (1:n_species) - 1]
          i_param <- i_param + n_species
          
          
        } else {
          
          # otherwise, without partial species spooling, 
          # fit beta independently between species 
          # Thus set mean to 0 and SD to 1 (i.e. latent variable z become directly the fitted parameter)
          p_mean_dbh <- 0
          p_dbh_per_sp <- p[i_param + (1:n_species) - 1]
          i_param <- i_param + n_species
        }
        
      } else {
        
        # single dbh effect over all species
        p_mean_dbh <- 0
        p_dbh_per_sp <- rep(p[i_param], n_species)
        i_param <- i_param + 1
      }
      
    } else {
      
      # Otherwise, no dbh effect
      p_mean_dbh <- 0
      p_dbh_per_sp <- rep(0, n_species)
    }
    
    
    # 4. Covariance between intercept/slope ----
    if (mod_design$consider_covariance) {
      
      # Rho parameter (correlation factor)
      p_rho_raw <- p[i_param]
      i_param <- i_param + 1
      
    } else {
      
      p_rho_raw <- 0
    }
    
    
    # List with all different group of parameters ----
    return(list(
      
      # Model SD residuals
      "sigma_log" = p_sigma,
      
      # Hierarchical random effect
      "sigma_origin_log" = p_sigma_origin,
      "sigma_site_log" = p_sigma_site,
      "z_random_per_origin" = p_z_origins,
      "z_random_per_site" = p_z_sites,
      
      # Predictors effect
      "intercept_per_sp" = p_intercept_per_sp,
      "dbh_per_sp" = p_dbh_per_sp,
      
      # Hyperparameters of predictors
      "mean_intercept" = p_mean_intercept,
      "mean_dbh" = p_mean_dbh,
      "sigma_intercept_log" = p_sigma_intercept,
      "sigma_dbh_log" = p_sigma_dbh,
      
      # Correlation parameter for MVN
      "rho_raw" = p_rho_raw
      
    ))
  }
  
  
  ## FUNCTION TO COMPUTE LOG LIKELIHOOD OF THE DATA ----
  compute_log_likelihood <- function(p, pointwise = FALSE, print.pb = FALSE) {
    
    
    # 1. Parameters of the LAD model ---- 
    
    ## 1.1. Unpack the parameter vector ----
    p_list <- split_parameters_vector(p)
  
    ## 1.2. Check for invalid parameter values ----
    # if (any(!is.finite(unlist(p_list)))) return(-Inf)

    
    
    # 2. Compute the hierarchical random effects ----
    # non-centered parameterization
    
    ### 2.1. Unlog sigma parameters ----
    p_sigma_origin <- exp(p_list$sigma_origin_log)
    p_sigma_site <- exp(p_list$sigma_site_log)
    
    ## 2.2. Compute the origins mean ----
    mean_origins <- p_sigma_origin * p_list$z_random_per_origin  # z_origin ~ N(0,1)
    
    ## 2.3. Get the associated origin mean for each site ----
    mean_origins_per_site <- mean_origins[id_origin_per_site]
    
    ## 2.4. Site effect nested in the origin effect ----
    p_random_per_site <- mean_origins_per_site + p_sigma_site * p_list$z_random_per_site  # z_site ~ N(0,1)
    
    
    
    # 3. Compute intercept and dbh effects ----
    
    ## 3.1. Transform back the correlation factor ----
    # Here, fit an unconstrained rho and apply a transformation to ensure -1 < rho < 1
    # ensure rho constrains with a better MCMC mixing and stability
    # plogis(x) in ]0,1[, 2*plogis(x) in ]0,2[ and 2*plogis(x) - 1 in ]-1,1[
    rho <- 2*plogis(p_list$rho_raw) - 1 
    
    
    ## 3.2. Transform back the SDs ----
    # To ensure SD positivity without constrained, by fitting in an infinite range
    sigma_intercept <- exp(p_list$sigma_dbh_log)
    sigma_dbh <- exp(p_list$sigma_dbh_log)
    
    
    ## 3.3. Compute the species parameters alpha and beta ----
    # Here, z_alpha/z_beta is a vector of latent variables for each species
    # Be careful, it is well ordered to respect the species order
    # z_alpha and z_beta ~ N(0,1)
    # If MVN, we compute SIGMA with the Cholesky factor L (SIGMA = L.t(L))
    p_intercept_per_sp <- p_list$mean_intercept + sigma_intercept * p_list$intercept_per_sp
    
    sqrt1mr2 <- sqrt(1 - rho^2) # Intermediate variables to faster the computation time
    p_dbh_per_sp <- p_list$mean_dbh + sigma_dbh * (rho * p_list$intercept_per_sp + sqrt1mr2 * p_list$dbh_per_sp)
    
    
    
    # 4. Compute the residuals with the given lad parameters ----
    
    ## 4.1. get residuals by running SmasaraLight ----
    out_sl <- compute_pacl_residuals(data_sensors,
                                     data_plots,
                                     data_calib, 
                                     data_rad, 
                                     p_random_per_site,
                                     p_intercept_per_sp,
                                     p_dbh_per_sp,
                                     # p_list$compet,
                                     print.pb)
    residuals <- out_sl$residuals
    
    ## 4.2. Check for invalid residuals
    # if (any(!is.finite(residuals))) return(-Inf)

    
    
    # 5. Compute the log-likelihood
    
    ## 5.1. Unlog sigma parameter ----
    p_sigma <- exp(p_list$sigma_log)
    
    
    ## 5.2. Total log-likelihood of the data ----
    # During the MCMC sampling,
    # We combine this log-likelihood with the log-priors computed by BayesianTools from priors defined in the BAyesian Setup 
    log_likelihood_data <- dnorm(residuals, 
                                 mean = 0, sd = p_sigma, 
                                 log = TRUE)
    
    ## 5.3. Check for invalid log-likelihood ----
    # if (any(!is.finite(log_likelihood_data))) return(-Inf)
  
    
    
    # 6. Return correct log-likelihood ----
    
    ## 6.1. Pointwise log-likelihood and residuals ----
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
    
    
    ## 6.2. Otherwise, return the log-likelihood scalar ----
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
  
  
  
  ## Standardized variables ----
  
  ### DBH
  dbh_vect <- data_calib %>% 
    purrr::map(~.x$trees$dbh_cm) %>% 
    purrr::reduce(c)
  
  dbh_mean <- mean(dbh_vect)
  dbh_sd <- sd(dbh_vect)
  
  data_calib <- data_calib %>% 
    purrr::map(~{
      .x$trees <- .x$trees %>% 
        dplyr::mutate(dbhstd = (dbh_cm - dbh_mean) / dbh_sd)
      .x
    })
  
  
  
  ## Global variables ----
  origin_names <- unique(data_plots$origin)
  site_names <- data_plots$name
  
  site_origins <- data_plots$origin # Origins of all the site
  id_origin_per_site <- match(site_origins, origin_names) # Corresponding id in origin_names for all the sites
  
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
  
  #### Intercept species partial pooling ----
  if (mod_design$intercept_sp_pooling) {
    par_normal <- c(par_normal, "sigma_intercept_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.5)
  }
  
  #### Dbh effect species partial pooling ----
  if (mod_design$dbh_sp_pooling) {
    par_normal <- c(par_normal, "sigma_dbh_log")
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
 
  
  ### Random effect latent variables ----
  # Non-centered parameterization with latent variables N(0,1)
  # To increase computation time and favour mixing
  
  #### Origin random effect ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, paste0("z_origin.", origin_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_origins))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_origins))
  }
  
  #### Site random effect ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, paste0("z_site.", site_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_sites))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_sites))
  }
  
  
  ### LAD intercept ----
  # (based on the best LAD guess from the methodology 1)
  # Here, we estimate the log(LAD), to ensure for LAD positivity
  lad_guess_method1 <- mean(output_lad_method1$best_lad[output_lad_method1$converged])
  
  if (lad_guess_method1 <= 0) stop("best lad guess from methodo is lower or equal to 0")
  lad_guess_method1_log <- log(lad_guess_method1)
  
  
  # Species-specific intercept
  if (mod_design$intercept_per_sp) {
    
    # Model the species-specific intercept within a normal distribution
    # centered around a mean species mean_intercept with a sd parameter sigma_intercept
    # the sigma paramter sigma_intercept has been defined previously
    if (mod_design$intercept_sp_pooling) {
      
      # Mean species with a prior centered around the predicted LAD from methodology 1
      par_normal <- c(par_normal, "mean_intercept")
      prior_normal_mean <- c(prior_normal_mean, lad_guess_method1_log)
      prior_normal_sd <- c(prior_normal_sd, 0.5) 
      
      # Latent variables (for species pooling)
      # To increase computation time and favour mixing: use a non-centered parametrization
      # Get the species intercept with alpha_s = mean + sigma * z_s (with z in N(0,1))
      par_normal <- c(par_normal, paste0("z_intercept.", species2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
      prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_species))
      
      
    } else {
      
      # fit intercept directly otherwise
      par_normal <- c(par_normal, paste0("intercept.", species2calib)) 
      prior_normal_mean <- c(prior_normal_mean, rep(lad_guess_method1_log, times=n_species))
      prior_normal_sd <- c(prior_normal_sd, rep(0.5, times=n_species))
    }
    
  } else {
    
    # a single intercept for all species
    par_normal <- c(par_normal, "intercept")
    prior_normal_mean <- c(prior_normal_mean, lad_guess_method1_log)
    prior_normal_sd <- c(prior_normal_sd, 0.5) 
  }
  
  
  
  ### DBH effect ----
  # BE CAREFUL : DBH in this model is standardized
  # LAD = intercept at DBHstd = 0 <=> DBH = mean(dataset)
  if (mod_design$dbh_effect) {
    
    if (mod_design$dbh_per_sp) {
      
      # Model the species-specific dbh effect within a normal distribution
      # centered around a single mean mean_dbh with a sd parameter sigma_dbh
      # the sigma paramter sigma_intercept has been defined previously
      if (mod_design$dbh_sp_pooling) {
        
        # Mean species dbh effect with a prior centered around 0 (no effect)
        par_normal <- c(par_normal, "mean_dbh")
        prior_normal_mean <- c(prior_normal_mean, 0)
        prior_normal_sd <- c(prior_normal_sd, 0.2) 
        
        # Latent variables (for species pooling)
        # To increase computation time and favour mixing: use a non-centered parametrization
        # Get the species dbh effect with alpha_s = mean + sigma * z_s (with z in N(0,1))
        par_normal <- c(par_normal, paste0("z_dbh.", species2calib))
        prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
        prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_species))
        
        
      } else {
        
        # fit dbh coefficient directly otherwise
        par_normal <- c(par_normal, paste0("dbh.", species2calib)) 
        prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
        prior_normal_sd <- c(prior_normal_sd, rep(0.2, times=n_species))
      }
      
    } else {
      
      # a single intercept for all species
      par_normal <- c(par_normal, "dbh")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.2) 
    }
    
  } # Otherwise, no dbh effect

  
  ### Covariance intercept/slope ----
  if (mod_design$consider_covariance) {
    
    # Send error if cannot model covariance 
    if (!mod_design$intercept_per_sp | 
        !mod_design$dbh_effect | 
        !mod_design$dbh_per_sp | 
        !mod_design$intercept_sp_pooling | 
        !mod_design$dbh_sp_pooling) {
      stop("To model the variance-covariance matrix, you need to consider a partial pooling species-specific intercept and dbh")
    }
    
    # Add rho parameter (correlation factor)
    # Here we model rho_raw which is back transformed into rho with rho = 2*plogis(rho_raw) - 1
    # Easier to mix for MCMC because unconstrained parameter, whereas -1 < rho < 1
    # And transformation transform rho from ]-Inf, +Inf[ to ]-1, 1[
    par_normal <- c(par_normal, "rho_raw")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 1) 
  }
  
  
  ## Competition effect
  if (mod_design$compet_effect) {
    stop("DEBUG: competition effect do not work for now")
    
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
    names = par_names,
    parallel = parallel
  )
  
  
  ## Return the environment ----
  return(rlang::current_env())
}
