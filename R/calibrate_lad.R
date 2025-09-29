calibrate_lad <- function(init_db, data_sl, data_rad, exp_design,
                          species2calib, site_names,
                          n_chains, n_iterations, n_burning,
                          output_folder) {
  
  # SUB-FUNCTIONS ----
  
  ## FUNCTION TO RUN SAMSARALIGHT ----
  run_sl_standXlad <- function(data_sl, 
                               data_rad, 
                               site, 
                               rep,
                               lad_sites,
                               lad_species,
                               lad_dbh) {
    
    # Create table of intercepts for each species
    # Here the order of the species are really important
    # It is the same as in the species2calib and in p_lad_intercepts
    sp_intercepts <- data.frame(
      species_calib = species2calib,
      sp_intercept = lad_species
    )
    
    # Get the random effect of the site
    # Here, the order is really important, and is the same between site_names and lad_sites vectors
    site_intercept <- lad_sites[which(site_names == site)]
    
    # Get tree dataset and set the LAD value for all trees
    tmp_trees <- data_sl[[site]][[rep]]$trees %>% 
      dplyr::left_join(sp_intercepts, by = "species_calib") %>% 
      dplyr::mutate( 
        crown_lad = site_intercept + sp_intercept + dbh_cm * lad_dbh)
    
    # Run SamsaraLight
    tmp_out_samsalight <- 
      sl_run(tmp_trees, 
             data_rad[[site]],
             sensors = data_sl[[site]][[rep]]$sensors, 
             sensors_only = TRUE,
             latitude = 46, slope = 0, 
             aspect = 144, north_to_x_cw = 54,
             start_day = 121, end_day = 273,
             cell_size = data_sl[[site]][[rep]]$info$cell_size, 
             n_cells_x = data_sl[[site]][[rep]]$info$n_cells_x, 
             n_cells_y = data_sl[[site]][[rep]]$info$n_cells_y,
             turbid_medium = TRUE,
             trunk_interception = FALSE,
             soc = TRUE,
             height_anglemin = 15,
             direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
             direct_anglestep = 5,
             diffuse_anglestep = 15)
    
    # Sl sensors output for the plot
    tmp_out_samsalight$sensors
    
  }
  
  ## FUNCTION TO RUN SAMSARALIGHT ON ALL STANDS ----
  compute_pacl_residuals <- function(data_sl,
                                     data_rad,
                                     data_sensors,
                                     exp_design,
                                     lad_sites,
                                     lad_species,
                                     lad_dbh) {
    
    # Output list
    out_residuals_list <- vector("list", nrow(exp_design))
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = nrow(exp_design),
                         style = 3, width = 50, char = "=")
    
    # Run for each simulation of the experimental design
    for (i in 1:nrow(exp_design)) {
      
      # Get info about the simulation
      tmp_simu_site <- as.character(exp_design$site[i])
      tmp_simu_rep <- exp_design$replicate[i]
      
      # Sl sensors output of all the plots
      tmp_out_sl <- run_sl_standXlad(data_sl, 
                                     data_rad, 
                                     tmp_simu_site, 
                                     tmp_simu_rep,
                                     lad_sites,
                                     lad_species,
                                     lad_dbh)
      
      # Compute mean residuals between all sensors
      out_residuals_list[[i]] <- dplyr::left_join(
        
        # Estimated pacl from virtual sensors
        tmp_out_sl %>%
          dplyr::select(id_sensor, pacl_sl = pacl_slope),
        
        # Measured pacl from field sensor
        data_sensors[[tmp_simu_site]] %>%
          dplyr::select(id_sensor = id, pacl_field = PACLtotal),
        
        by = "id_sensor"
      ) %>% 
        
        dplyr::mutate(residuals = pacl_sl - pacl_field) %>% 
        dplyr::select(id_sensor, pacl_sl, pacl_field, residuals)
      
      
      # Update the progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    dplyr::bind_rows(out_residuals_list)
  }
  
  
  
  ## FUNCTION TO COMPUTE LOG POSTERIOR FROM LIKELIHOODS ----
  compute_log_posterior <- function(p) {
    
    # 1. Parameters for the calibration run
    p_sigma <- p[1] # Standard deviation of residuals
    p_sigma_sites <- p[2] # SD of the site effect
    p_lad_sites <- p[(1:n_sites) + 2] #  Site intercepts
    p_lad_species <- p[(1:n_species) + 2 + n_sites] # Species intercepts
    p_lad_dbh <- p[1 + (2 + n_sites + n_species)] # DBH effect 
    
    
    # 2. Compute the residuals of the given lad parameters
    out_sl <- compute_pacl_residuals(data_sl, 
                                     data_rad, 
                                     init_db$sensors,
                                     exp_design, 
                                     p_lad_sites,
                                     p_lad_species,
                                     p_lad_dbh)
    residuals <- out_sl$residuals
    
    # 3. Compute the log-likelihood of the hirerarchcal random site effect
    log_likelihood_sites <- sum(dnorm(p_lad_sites, mean = 0, sd = p_sigma_sites, log = TRUE)) 
    
    # 4. Compute the data log-likelihood
    log_likelihood_data <- sum(dnorm(residuals, mean = 0, sd = p_sigma, log = TRUE)) 
    
    # 5. Copute the total log likelihood
    log_posterior <- log_likelihood_sites + log_likelihood_data
    
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
    best <- c(halfcauchy_scale, (uniform_upper - uniform_lower) / 2)
    
    return(BayesianTools::createPrior(
      density = prior_density,
      sampler = prior_sampler,
      lower = lower,
      upper = upper,
      best = best
    ))
  }
  
  
  # SCRIPT ----
  
  ## Base variables ----
  n_species <- length(species2calib)
  n_sites <- length(site_names)

  
  ## Prior definition ----
  
  # HALFCAUCHY PRIORS MUST BE THE FIRST ONES (typically for modelling SD)
  par_names <- c(
    "sigma", # model sigma SD
    "sigma_site", # SD of the random site effect (hyperparameter)
    paste0("site.", site_names), # random sites intercept
    paste0("species.", species2calib), # species intercept
    "dbh" # dbh effect
  )
  
  halfcauchy_S <- rep(5, times=2)
  
  uniform_LB <- c(rep(-10, times=n_sites), 
                  rep(0.1, times=n_species), 
                  -0.001)
  
  uniform_UB <- c(rep(10, times=n_sites),
                  rep(2, times=n_species),
                  0.001) 
  
  prior <- createCombinedPriors(halfcauchy_scale = halfcauchy_S,
                                uniform_lower = uniform_LB, 
                                uniform_upper = uniform_UB) 
  
  
  
  ## Bayesian setup ----
  bayesianSetup <- BayesianTools::createBayesianSetup(compute_log_posterior, 
                                                      prior, 
                                                      names = par_names)
  
  settings <- list(iterations = n_iterations, nCR = n_chains, gamma = NULL, eps = 0, e = 0.05, 
                   pCRupdate = FALSE, updateInterval = 100, burnin = n_burning, thin = 1, 
                   adaptation = 0.2, parallel = NULL, Z = NULL, ZupdateFrequency = 10, 
                   pSnooker = 0.1, DEpairs = 2, consoleUpdates = 100, startValue = NULL, 
                   currentChain = 1, message = TRUE)
  
  
  ## Run Bayesian calibration ----
  t_start <- Sys.time()
  out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup, 
                                sampler = "DREAMzs", 
                                settings = settings)
  t_end <- Sys.time()
  
  computation_time <- t_end - t_start
  computation_time
  
  
  ## Save Bayesian output and configuration ----
  session <- Sys.info()
  thetime <- format(Sys.time(), "%Y_%m_%d_%Hh%M")
  output_fp <- file.path(output_folder, paste0("res_",thetime,".Rdata"))
  
  save(out,prior,bayesianSetup,settings,session,computation_time,
       file = output_fp)
  
  ## Return the output ----
  return(
    list(
      "output" = out,
      "info_calib_fp" = output_fp
    )
  )
}
