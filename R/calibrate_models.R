calibrate_models <- function(models_setup,
                             sampling_algo,
                             id_model,
                             i_chain,
                             logs_folder) {
  
  # Function to calibrate a model given the environment is set up ----
  calibrate_model <- function(sampling_algo) {
    
    # Define settings ----
    settings <- list(iterations = mod_design$n_iterations, 
                     nrChains = 1, burnin = 0, 
                     gamma = NULL, eps = 0, e = 0.05, 
                     pCRupdate = TRUE, updateInterval = 50, thin = 1, 
                     adaptation = 0.2, Z = NULL, ZupdateFrequency = 10, 
                     pSnooker = 0.3, DEpairs = 2, startValue = NULL, 
                     consoleUpdates = 1, message = TRUE)
    
    
    
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
  
  
  
  # Initialize the logs ----
  
  ## Create the folder
  dir.create(logs_folder, showWarnings = FALSE)
  
  # Unique log file per worker that close at the end of the function
  log_file <- sprintf("logs/calib-mod%s-chain%i-worker%s-%s.log", 
                      id_model, i_chain, 
                      Sys.getpid(), Sys.Date())
  
  # Capture output (prints)
  con <- file(log_file, open = "wt")
  sink(con)

  
  # Calibrate the model ----
  ## Run the model by storing errors, warnings and messages
  cat("----- Calibrating model", id_model, "/ chain", i_chain, "-----\n\n\n")
  
  result <- tryCatch(
    {
      # Use withCallingHandlers to capture messages/warnings
      withCallingHandlers({
        
        cat("Starting MCMC\n\n")
        
        ## Set the model environment
        environment(calibrate_model) <- models_setup[[id_model]]
        
        ## Calibrate the model
        out <- calibrate_model(sampling_algo)
        
        ## Add model IDs to output
        out$id_model <- id_model
        out$i_chain <- i_chain
        
        cat("\n\nFinished MCMC\n\n")
        
        out
      },
      message = function(m) {
        cat("\n[MESSAGE] ", conditionMessage(m), "\n", file = con)
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        cat("\n[WARNING] ", conditionMessage(w), "\n", file = con)
        invokeRestart("muffleWarning")
      })
    },
    error = function(e) {
      cat("\n[ERROR] ", conditionMessage(e), "\n", file = con)
      # Return a safe placeholder (e.g., NA or empty list)
      list(failed = TRUE, params = params, result = NA)
    }
  )
  
  cat("\n\n----- Finished model", id_model, "/ chain", i_chain, "-----")
  close(con)
  
  
  # Return the list of models' output ----
  return(result) 
}
