calibrate_models <- function(models_setup,
                             id_model,
                             i_rep,
                             logs_folder) {
  
  # Function to calibrate a model given the environment is set up ----
  calibrate_model <- function() {
    
    
    # Define settings ----
    settings <- list(iterations = mod_design$n_iterations,
                     burnin = mod_design$n_burning,
                     thin = 1,
                     nCR = mod_design$n_subchains,
                     nrChains = mod_design$n_chains, 
                     gamma = NULL, 
                     eps = 1e-6, e = 0.05, 
                     pCRupdate = TRUE, updateInterval = 50, 
                     adaptation = 0.2, 
                     Z = NULL, ZupdateFrequency = 10, 
                     pSnooker = 0.1, DEpairs = 3, 
                     startValue = NULL, 
                     consoleUpdates = 1, message = TRUE)
    
    
    
    # Run Bayesian calibration ----
    t_start <- Sys.time()
    out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup, 
                                  sampler = "DREAMzs", 
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
  
  # Create the folder
  dir.create(logs_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Unique log file per worker that close at the end of the function
  log_file <- sprintf(file.path(logs_folder,"calib-mod%s-chain%i-worker%s-%s.log"), 
                      id_model, i_rep, 
                      Sys.getpid(), Sys.Date())
  
  
  # Redirect ALL stdout + C output
  log_con <- file(log_file, open = "at")
  
  sink(log_con, type = "output", append = TRUE)
  sink(log_con, type = "message", append = TRUE)
  
  # Close connections at the end of the function
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
  }, add = TRUE)
  
  
  ### Manual logger
  # R does not immediately write to disk.
  # It puts the text in a buffer in memory
  # Thus, flush() pushes everything in the buffer to disk
  log_write <- function(...) {
    cat(..., "\n")
    flush.console()
  }
  
  
  # Calibrate the model ----
  ## Run the model by storing errors, warnings and messages
  log_write("----- Calibrating model", id_model, "/ rep", i_rep, "-----\n")
  
  result <- tryCatch(
    {
      # Use withCallingHandlers to capture messages/warnings
      withCallingHandlers({
        
        log_write("Starting MCMC\n")
        
        ## Set the model environment
        environment(calibrate_model) <- models_setup[[id_model]]
        
        ### capture printed output
        out <- calibrate_model()
        
        ## Add model IDs to output
        out$id_model <- id_model
        out$i_rep <- i_rep
        
        log_write("\nFinished MCMC\n")
        
        out
      },
      message = function(m) {
        log_write("[MESSAGE]", conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        log_write("[WARNING]", conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    },
    error = function(e) {
      log_write("[ERROR]", conditionMessage(e))
      # Return a safe placeholder (e.g., NA or empty list)
      list(failed = TRUE, error = conditionMessage(e))
    }
  )
  
  log_write("\n----- Finished model", id_model, "/ chain", i_rep, "-----\n")
  
  
  # Return the list of models' output ----
  return(result) 
}