initialise_models <- function(exp_design,
                              data_plots,
                              data_stands,
                              data_rad,
                              output_lad_method1,
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
                                                 output_lad_method1,
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
                             output_lad_method1,
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
                               lad_dbh_per_sp,
                               lad_batot_per_sp,
                               lad_dbhXbatot_per_sp
                               ) {
    
    # Species-specific parameters ----
    # Here the order of the species modalities is really important
    # It is the same as in the species2calib, in lad_intercept_per_species...
    # And also the same in is_species_gymno when we added gymnosperm coefficients to species-specific parameters
    sp_coefs <- data.frame(
      species = species2calib,
      intercept = lad_intercept_per_sp,
      beta_dbh = lad_dbh_per_sp,
      beta_batot = lad_batot_per_sp,
      beta_dbhxbatot = lad_dbhXbatot_per_sp
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
          beta_dbh * dbhstd + 
          beta_batot * batotstd + 
          beta_dbhxbatot * dbhstd * batotstd,
        crown_lad = exp(eta)
      )
    
    
    # Run SamsaraLight ----
    tmp_out_samsalight <- SamsaRaLight::run_sl(
      tmp_stand, 
      data_rad[[site]],
      sensors_only = TRUE,
      turbid_medium = TRUE,
      use_torus = TRUE,
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
                                     lad_dbh_per_sp,
                                     lad_batot_per_sp,
                                     lad_dbhXbatot_per_sp,
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
        lad_dbh_per_sp,
        lad_batot_per_sp,
        lad_dbhXbatot_per_sp
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
    
    
    # 2. Site/origin hierarchical random effects ----
    
    ## 2.1. Origin latent variables ----
    if (mod_design$origin_rd_effect) {
      p_z_origins <- p[i_param + (1:n_origins) - 1]
      i_param <- i_param + n_origins
    } else {
      p_z_origins <- rep(0, n_origins)
    }
    
    ## 2.2. Site latent variables ----
    if (mod_design$site_rd_effect) {
      p_z_sites <- p[i_param + (1:n_sites) - 1]
      i_param <- i_param + n_sites
    } else {
      p_z_sites <- rep(0, n_sites)
    }
  
    
    
    # 3. Predictors effect ----
    
    ## 3.1. Intercept ----
    
    # Mean species intercept
    p_intercept <- p[i_param]
    i_param <- i_param + 1 
    
    # Species pooling latent variable for intercept
    if (mod_design$intercept_sp_pooling) {
      
      p_z_intercept_sp <- p[i_param + (1:n_species) - 1]
      i_param <- i_param + n_species
    } else {
      p_z_intercept_sp <- rep(0, n_species)
    }
    
    # Intercept for gymnosperms
    if (mod_design$group_effect) {
      
      p_intercept_gymno <- p[i_param]
      i_param <- i_param + 1
    } else {
      p_intercept_gymno <- 0 
    }
    
    
    
    ## 3.2. DBH effect ----
    if (mod_design$dbh_effect) {
      
      # Mean species dbh
      p_dbh <- p[i_param]
      i_param <- i_param + 1 
      
      # Species pooling latent variables for dbh
      if (mod_design$dbh_sp_pooling) {
        
        p_z_dbh_sp <- p[i_param + (1:n_species) - 1]
        i_param <- i_param + n_species
      } else {
        p_z_dbh_sp <- rep(0, n_species)
      }
      
      # DBH effect for gymnosperms
      if (mod_design$dbh_per_group) {
        
        p_dbh_gymno <- p[i_param]
        i_param <- i_param + 1
      } else {
        p_dbh_gymno <- 0 
      }
      
    } else {
      
      # Otherwise, no dbh effect
      p_dbh <- 0
      p_z_dbh_sp <- rep(0, n_species) 
      p_dbh_gymno <- rep(0, n_species)
    }
    
    
    
    ## 3.3. Covariance between intercept/slope ----
    if (mod_design$covariance_sp_pooling) {
      
      # Rho parameter (correlation factor)
      p_rho_raw <- p[i_param]
      i_param <- i_param + 1
      
    } else {
      
      p_rho_raw <- 0
    }
    
    
    ## 3.4. Within-site effect of total basal area ----
    if (mod_design$batot_site_effect) {
      
      # Mean species batot effect
      p_batot <- p[i_param]
      i_param <- i_param + 1 
      
      # BATOT effect for gymnosperms
      if (mod_design$batot_per_group) {
        
        p_batot_gymno <- p[i_param]
        i_param <- i_param + 1
      } else {
        p_batot_gymno <- 0 
      }
      
    } else {
      
      # Otherwise, no dbh effect
      p_batot <- 0
      p_batot_gymno <- rep(0, n_species)
    }
    
    
    
    ## 3.5. BATOTxDBH interaction ----
    if (mod_design$dbh_interaction_batot) {
      
      # Mean species batot effect
      p_dbhXbatot <- p[i_param]
      i_param <- i_param + 1 
      
      # BATOT effect for gymnosperms
      if (mod_design$interaction_per_group) {
        
        p_dbhXbatot_gymno <- p[i_param]
        i_param <- i_param + 1
      } else {
        p_dbhXbatot_gymno <- 0 
      }
      
    } else {
      
      # Otherwise, no dbh effect
      p_dbhXbatot <- 0
      p_dbhXbatot_gymno <- rep(0, n_species)
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
      
      # Intercept
      "intercept" = p_intercept,
      "z_intercept_species" = p_z_intercept_sp,
      "sigma_intercept_species_log" = p_sigma_intercept,
      "intercept_gymno" = p_intercept_gymno,
      
      # DBH
      "dbh" = p_dbh,
      "z_dbh_species" = p_z_dbh_sp,
      "sigma_dbh_species_log" = p_sigma_dbh,
      "dbh_gymno" = p_dbh_gymno,
      
      # Correlation parameter for MVN
      "rho_raw" = p_rho_raw,
      
      # BATOT
      "batot" = p_batot,
      "batot_gymno" = p_batot_gymno,
      
      # Interaction BATOTxDBH
      "dbhXbatot" = p_dbhXbatot,
      "dbhXbatot_gymno" = p_dbhXbatot_gymno
      
      
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
    
    
    
    # 2. Compute the hierarchical random effects ----
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
    
    
    # 3. Compute species-pooling intercept and dbh effects with MVN ----
    
    ## 3.1. Transform back the correlation factor ----
    # Here, fit an unconstrained rho and apply a transformation to ensure -1 < rho < 1
    # ensure rho constrains with a better MCMC mixing and stability
    # plogis(x) in ]0,1[, 2*plogis(x) in ]0,2[ and 2*plogis(x) - 1 in ]-1,1[
    rho <- 2*plogis(p_list$rho_raw) - 1 
    
    
    ## 3.2. Transform back the SDs ----
    # To ensure SD positivity without constrained, by fitting in an infinite range
    sigma_intercept <- exp(p_list$sigma_intercept_species_log)
    sigma_dbh <- exp(p_list$sigma_dbh_species_log)
    
    
    ## 3.3. Compute the species parameters alpha and beta ----
    # Here, z_alpha/z_beta is a vector of latent variables for each group modalities
    # Be careful, it is well ordered to respect the group modalities order
    # z_alpha and z_beta ~ N(0,1)
    # If MVN, we compute SIGMA with the Cholesky factor L (SIGMA = L.t(L))
    p_intercept_per_sp <- p_list$intercept + sigma_intercept * p_list$z_intercept_species
    
    sqrt1mr2 <- sqrt(1 - rho^2) # Intermediate variables to faster the computation time
    p_dbh_per_sp <- p_list$dbh + sigma_dbh * (rho * p_list$z_intercept_species + sqrt1mr2 * p_list$z_dbh_species)
    
    
    # 4. Add the effect of functional group for each species ----
    # Add the effect only for gymnoseprm species (angiosperm is the reference group)
    p_intercept_per_sp <- p_intercept_per_sp + is_species_gymno * p_list$intercept_gymno
    p_dbh_per_sp <- p_dbh_per_sp + is_species_gymno * p_list$dbh_gymno
    p_batot_per_sp <- rep(p_list$batot, times = n_species) + is_species_gymno * p_list$batot_gymno
    p_dbhXbatot_per_sp <- rep(p_list$dbhXbatot, times = n_species) + is_species_gymno * p_list$dbhXbatot_gymno
    
    
    # 5. Compute the residuals with the given lad parameters ----
    # (difference between estimated/observed PACL for each sensor)
    
    ## 5.1. get residuals by running SamsaraLight ----
    out_sl <- compute_pacl_residuals(data_stands, 
                                     data_rad, 
                                     p_random_per_site,
                                     p_intercept_per_sp,
                                     p_dbh_per_sp,
                                     p_batot_per_sp,
                                     p_dbhXbatot_per_sp,
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
  # Remove sensors that did not converged
  
  if (mod_design$filter_sensors) {
    
    for (site_name in names(data_stands)) {
      
      data_stands[[site_name]]$sensors <- data_stands[[site_name]]$sensors %>%
        dplyr::left_join(
          output_lad_method1 %>% 
            dplyr::filter(site == site_name) %>%
            dplyr::ungroup() %>% 
            dplyr::select(id_sensor, converged),
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
  
  ### DBH (per trees) ----
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
  
  
  ### BATOT (per site) ----
  batot_vect <- data_stands %>% 
    purrr::map(~unique(.x$trees$batot_m2ha)) %>% 
    purrr::reduce(c)
  
  batot_mean <- mean(batot_vect)
  batot_sd <- sd(batot_vect)
  
  data_stands <- data_stands %>% 
    purrr::map(~{
      .x$trees <- .x$trees %>% 
        dplyr::mutate(batotstd = (batot_m2ha - batot_mean) / batot_sd)
      .x
    })
  
  
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
  
  
  ### Fonctional group effect ----
  group_per_species <- data_stands %>%
    purrr::map(~.x$trees %>% 
                 dplyr::select(species, functional_group) %>% 
                 dplyr::distinct()) %>% 
    dplyr::bind_rows() %>% 
    dplyr::distinct() %>% 
    # Order species with the same order as the species vector
    dplyr::mutate(species = factor(species, levels = species2calib)) %>% 
    dplyr::arrange() %>% 
    dplyr::pull(functional_group)
  
  is_species_gymno <- group_per_species == "gymnosperm"
  
  
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
  prior_normal_mean <- c(prior_normal_mean, log(0.05))
  prior_normal_sd <- c(prior_normal_sd, 0.1)
  
  
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
  
  
  #### Intercept species partial pooling ----
  if (mod_design$intercept_sp_pooling) {
    par_normal <- c(par_normal, "sigma_intercept_sp_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.3)
  }
  
  #### Dbh effect species partial pooling ----
  if (mod_design$dbh_sp_pooling) {
    par_normal <- c(par_normal, "sigma_dbh_sp_log")
    prior_normal_mean <- c(prior_normal_mean, log(0.3))
    prior_normal_sd <- c(prior_normal_sd, 0.3)
  }
  
  
  ### Hierarchical site/origin random effect ----
  # Non-centered parameterization with latent variables N(0,1)
  # To increase computation time and favour mixing
  
  #### Origin latent variables ----
  if (mod_design$origin_rd_effect) {
    par_normal <- c(par_normal, paste0("z_origin.", origin_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_origins))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_origins))
  }
  
  #### Site latent variables ----
  if (mod_design$site_rd_effect) {
    par_normal <- c(par_normal, paste0("z_site.", site_names))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_sites))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_sites))
  }
  
  
  ### Mean species LAD intercept ----
  # Here, we estimate the log(LAD), to ensure for LAD positivity

  if (prior_lad <= 0) stop("prior lad is lower or equal to 0")
  prior_lad_log <- log(prior_lad)
  

  #### Mean species intercept ----
  par_normal <- c(par_normal, "intercept")
  prior_normal_mean <- c(prior_normal_mean, prior_lad_log)
  prior_normal_sd <- c(prior_normal_sd, 0.5)
  
  
  #### Species-pooling ----
  # Model the species-specific intercept within a normal distribution
  # centered around the mean species intercept intercept with a sd parameter sigma_intercept
  # the sigma paramter sigma_intercept has been defined previously
  # Same as modelling an intercept with species as random effect
  if (mod_design$intercept_sp_pooling) {
    
    # Latent variables (for species pooling)
    # To increase computation time and favour mixing: use a non-centered parametrization
    # Get the group intercept with alpha_s = mean + sigma * z_s (with z in N(0,1))
    par_normal <- c(par_normal, paste0("intercept_z_sp.", species2calib))
    prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_species))
    prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_species))
    
  } 
  
  
  #### Functional group intercept effect ----  
  # Intercept change compared to the ref group (ref is angiosperm)
  # Normal distribution centered around 0
  if (mod_design$group_effect) {
    
    par_normal <- c(par_normal, "intercept_gymno")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
    
  } # Otherwise, no group effect
  
  
  
  ### DBH effect ----
  # BE CAREFUL : DBH in this model is standardized
  # LAD = intercept at DBHstd = 0 <=> DBH = mean(dataset)
  if (mod_design$dbh_effect) {
    
    # Cannot have dbh species pooling with group-specific slope on dbh
    if (mod_design$dbh_per_group & mod_design$dbh_sp_pooling) {
      stop("DEBUG: cannot have dbh species pooling with group-specific slope on dbh")
    }
    
    ### Single slope for all species ----
    par_normal <- c(par_normal, "dbh")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
    
    
    #### Species pooling ----
    # Model the species-specific dbh effect within a normal distribution
    # centered around a mean dbh with a sd parameter sigma_dbh_sp
    # the sigma paramter sigma_dbh_sp has been defined previously
    if (mod_design$dbh_sp_pooling) {
      
      # Latent variables (for species pooling)
      # To increase computation time and favour mixing: use a non-centered parametrization
      # Get the group dbh effect with alpha_s = mean + sigma * z_s (with z in N(0,1))
      par_normal <- c(par_normal, paste0("dbh_z_sp.", groups2calib))
      prior_normal_mean <- c(prior_normal_mean, rep(0, times=n_groups))
      prior_normal_sd <- c(prior_normal_sd, rep(1, times=n_groups))
      
    } 
    
    ### Group-specific slope ----
    # Slope change compared to the ref group (ref is angiosperm)
    # Normal distribution centered around 0
    if (mod_design$dbh_per_group) {
      
      par_normal <- c(par_normal, "dbh_gymno")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
      
    }
    
  } # Otherwise, no dbh effect
  
  
  ### Covariance intercept/dbh when species pooling ----
  if (mod_design$covariance_sp_pooling) {
    
    # Send error if cannot model covariance 
    if (!mod_design$intercept_sp_pooling | 
        !mod_design$dbh_sp_pooling) {
      stop("DEBUG: to model the variance-covariance matrix, you need to consider a partial pooling species-specific intercept and dbh")
    }
    
    # Add rho parameter (correlation factor)
    # Here we model rho_raw which is back transformed into rho with rho = 2*plogis(rho_raw) - 1
    # Easier to mix for MCMC because unconstrained parameter, whereas -1 < rho < 1
    # And transformation transform rho from ]-Inf, +Inf[ to ]-1, 1[
    par_normal <- c(par_normal, "rho_raw")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 1) 
  }
  
  
  
  ### BATOT effect ----
  if (mod_design$batot_site_effect) {
    
    #### Single slope for all species ----
    par_normal <- c(par_normal, "batot")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
    
    #### Group-specific slope ----
    # (ref is angiosperm)
    if (mod_design$batot_per_group) {
      
      par_normal <- c(par_normal, "batot_gymno")
      prior_normal_mean <- c(prior_normal_mean, 0)
      prior_normal_sd <- c(prior_normal_sd, 0.5)
      
    }
    
  } # Otherwise, no batot effect
  
  
  ### Interaction BATOT x DBH ----
  if (mod_design$dbh_interaction_batot) {
    
    # Send error if interaction is set to true with dbh species pooling also activated
    if (mod_design$dbh_sp_pooling) {
      stop("DEBUG: DBH:BATOT does not work with dbh species pooling")
    }
    
    # Send error if batot or dbh effect are not activated
    if (!mod_design$dbh_effect |
        !mod_design$batot_site_effect) {
      stop("DEBUG: include both batot and dbh effect to consider interaction DBH:BATOT")
    }
    
    # Send error if batot or dbh effect are not considered equally (global or per group)
    if (mod_design$dbh_per_group != mod_design$batot_per_group) {
      stop("DEBUG: dbh and batot per group must be considered the same manner to consider batot:dbh interaction")
    }
    
    #### Single slope for all species ----
    par_normal <- c(par_normal, "dbhXbatot")
    prior_normal_mean <- c(prior_normal_mean, 0)
    prior_normal_sd <- c(prior_normal_sd, 0.5)
    
    #### Group-specific slope ----
    # (ref is angiosperm)
    if (mod_design$interaction_per_group) {
      
      par_normal <- c(par_normal, "dbhXbatot_gymno")
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
