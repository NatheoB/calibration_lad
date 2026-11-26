get_summary_pointwise_models <- function(models_setup,
                                         model_output,
                                         n_burning_per_rep,
                                         n_samples_per_rep,
                                         n_site_workers,
                                         logs_folder) {
  
  get_summary_pointwise_model <- function(output, n_burning_per_rep, n_samples_per_rep, n_site_workers) {
    
    # Start the cluster to run the sites in parallel ----
    if (n_site_workers > 1) {
      cl_sites <<- start_cluster_sites(parent.env(environment()), n_site_workers)
      on.exit({
        parallel::stopCluster(cl_sites)
        cl_sites <<- NULL
      }, add = TRUE)
    }
    
    # The run must have been calibrated with this setup (same parameters, same order of site effects)
    if (!identical(output$par_names, bayesianSetup$names) | !identical(output$site_names, site_names)) {
      stop("The run output does not match the model setup (parameter names or site order differ): recalibrate or rebuild the setup")
    }
    
    # Get samples of the internal chains of the run, after burn-in, with their chain identity
    chain_samples <- get_chain_samples(output, n_burning_per_rep, n_samples_per_rep)
    samples <- chain_samples$samples
    
    # Get the number of observation (i.e. total number of sesnors over all the sites)
    n_sensors <- nrow(data_stands %>% 
                        purrr::map(~.x$sensors) %>% 
                        dplyr::bind_rows())

    # Get the matrix of PACL residuals and pointwise log-likelihood for all samples of parameters
    out_residuals <- matrix(NA, nrow = nrow(samples), ncol = n_sensors)
    out_llpointwise <- matrix(NA, nrow = nrow(samples), ncol = n_sensors)
    
    pb <- txtProgressBar(min = 0, max = nrow(samples),
                         style = 3, width = 50, char = "=")

    for (i in 1:nrow(samples)) {

      # The important point here is the pointwise set to TRUE
      matrices_list <- compute_log_likelihood(samples[i,], pointwise = TRUE, print.pb = FALSE)
      
      out_residuals[i,] <- matrices_list$residuals
      out_llpointwise[i,] <- matrices_list$loglikelihood
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    list(
      "residuals" = out_residuals,
      "llpointwise" = out_llpointwise,
      "sensors" = matrices_list$sensors,
      "sigma" = exp(samples[, "sigma_mod_log"]),
      "id_chain" = chain_samples$id_chain
    )
  }
  
  # Get model infos ----
  id_model <- model_output$id_model
  i_rep <- model_output$i_rep
  
  
  
  # Initialize the logs ----
  
  # Create the folder
  dir.create(logs_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Unique log file per worker that close at the end of the function
  log_file <- sprintf(file.path(logs_folder,"pointwise-mod%s-chain%i-worker%s-%s.log"), 
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
  
  
  # Compute the residuals of all the last samples of each model ----
  ## Storing errors, warnings and messages
  log_write("----- Computing residuals and pointwise log-likelihood of model", id_model, "/ chain", i_rep, "-----\n\n\n")
  
  result <- tryCatch(
    {
      # Use withCallingHandlers to capture messages/warnings
      withCallingHandlers({
        
        log_write("Starting computing\n\n")
      
        ## Set the model environment
        environment(get_summary_pointwise_model) <- models_setup[[as.character(id_model)]]
        
        ## Get the model residuals and pointwise log-likelihood matrices
        out <- get_summary_pointwise_model(model_output, n_burning_per_rep, n_samples_per_rep, n_site_workers)
        
        ## Add model IDs to output
        out$id_model <- id_model
        out$i_rep <- i_rep
        
        log_write("\n\nFinished computing\n\n")
        
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
      list(failed = TRUE, result = NA)
    }
  )
  
  log_write("\n\n----- Finished model", id_model, "/ chain", i_rep, "-----")
  
  
  # Return the list of models' summary matrices ----
  return(result) 

  
}