calibrate_models <- function(models_setup,
                             n_chains,
                             n_iterations,
                             n_burning,
                             sampling_algo) {
  
  # Function to calibrate a model given the environment is set up ----
  calibrate_model <- function(settings, sampling_algo) {
    
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
  
  
  # Define settings ----
  settings <- list(iterations = n_iterations, burnin = n_burning, nrChains = n_chains,
                   gamma = NULL, eps = 0, e = 0.05, 
                   pCRupdate = TRUE, updateInterval = 50, thin = 1, 
                   adaptation = 0.2, parallel = NULL, Z = NULL, ZupdateFrequency = 10, 
                   pSnooker = 0.3, DEpairs = 2, consoleUpdates = 100, startValue = NULL, 
                   message = TRUE)
  
  
  # Calibrate all the models ----
  n_mods <- length(models_setup)
  models_output <- vector("list", n_mods)
  for (i in 1:n_mods) {
    
    # Print a message to follow the calibration process
    message(paste0("Calibrating model ", i, "/", n_mods, "..."))
    
    # Set the model environment
    environment(calibrate_model) <- models_setup[[i]]
    
    # Run the model
    models_output[[i]] <- calibrate_model(settings, sampling_algo)
  }
  
  
  # Return the list of models' output ----
  return(models_output)
}