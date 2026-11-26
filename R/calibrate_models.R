calibrate_models <- function(models_setup,
                             id_model,
                             i_rep,
                             n_site_workers,
                             logs_folder) {
  
  # Function to calibrate a model given the environment is set up ----
  calibrate_model <- function(n_site_workers) {
    
    # Start the cluster to run the sites in parallel ----
    if (n_site_workers > 1) {
      cl_sites <<- start_cluster_sites(parent.env(environment()), n_site_workers)
      on.exit({
        parallel::stopCluster(cl_sites)
        cl_sites <<- NULL
      }, add = TRUE)
    }
    
    
    # Define settings ----
    settings <- list(iterations = mod_design$n_iterations,
                     burnin = mod_design$n_burning,
                     thin = 1,
                     nCR = mod_design$nCR,
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
        "chains" = out$chain, # coda::mcmc.list, one element per internal DREAMzs chain (parameters + Lposterior, Llikelihood, Lprior); the bayesianSetup is not kept (its likelihood carries the whole setup environment, ~190 MB)
        "par_names" = bayesianSetup$names,
        "site_names" = site_names, # Order of the site random effects in the parameter vector
        "settings" = settings,
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
        environment(calibrate_model) <- models_setup[[as.character(id_model)]]
        
        ### capture printed output
        out <- calibrate_model(n_site_workers)
        
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