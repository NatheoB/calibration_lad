calibrate_models <- function(models_setup,
                             n_chains,
                             n_iterations,
                             n_burning,
                             sampling_algo) {
  
  # Function to calibrate a model given the environment is set up ----
  calibrate_model <- function(sampling_algo) {
    
    # Define settings ----
    settings <- list(iterations = mod_design$n_iterations, 
                     nrChains = mod_design$n_chains,
                     burnin = 0, 
                     gamma = NULL, eps = 0, e = 0.05, 
                     pCRupdate = TRUE, updateInterval = 50, thin = 1, 
                     adaptation = 0.2, parallel = NULL, Z = NULL, ZupdateFrequency = 10, 
                     pSnooker = 0.3, DEpairs = 2, consoleUpdates = 1, startValue = NULL, 
                     message = TRUE)
    
    # Run Bayesian calibration ----
    t_start <- Sys.time()
    out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup, 
                                  sampler = sampling_algo, 
                                  settings = settings)
    t_end <- Sys.time()
    
    # Set the computation time ----
    computation_time <- t_end - t_start
  
    # Return the output ----
    return(
      list(
        "outputs" = out,
        "computation_time" = computation_time
      )
    )
    
  }
  
  
  # Calibrate all the models ----
  ids_simu <- names(models_setup)
  n_mods <- length(ids_simu)
  
  models_output <- setNames(vector("list", n_mods), ids_simu)
  for (i in 1:n_mods) {
    
    # Get the simu name
    id_simu <- ids_simu[i]
    
    # Print a message to follow the calibration process
    message(paste0("Calibrating model ", id_simu, " - ", i, "/", n_mods, "..."))
    
    # Set the model environment
    environment(calibrate_model) <- models_setup[[id_simu]]
    
    # Run the model
    models_output[[id_simu]] <- calibrate_model(sampling_algo)
  }
  
  
  # Return the list of models' output ----
  return(models_output)
}
