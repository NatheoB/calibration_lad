get_summary_pointwise_models <- function(models_setup,
                                         model_output,
                                         logs_folder) {
  
  get_summary_pointwise_model <- function(output) {
    
    # Get samples of the chains of the model (bind the chains samplings with coda argument)
    samples <- BayesianTools::getSample(output, coda = F)
    samples <- samples[(nrow(samples)-mod_design$n_analysis+1):nrow(samples),]
    
    # Get the number of observation (i.e. total number of sesnors over all the sites)
    n_sensors <- nrow(data_sensors %>% dplyr::bind_rows())

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
      "llpointwise" = out_llpointwise
    )
  }
  
  # Get model infos ----
  id_model <- model_output$id_model
  i_chain <- model_output$i_chain
  
  
  # Initialize the logs ----
  
  ## Create the folder
  dir.create(logs_folder, showWarnings = FALSE)
  
  # Unique log file per worker that close at the end of the function
  log_file <- sprintf("logs/pointwise-mod%s-chain%i-worker%s-%s.log", 
                      id_model, i_chain, 
                      Sys.getpid(), Sys.Date())
  
  # Capture output (prints)
  con <- file(log_file, open = "wt")
  sink(con)
  
  
  # Compute the residuals of all the last samples of each model ----
  ## Storing errors, warnings and messages
  cat("----- Computing residuals and pointwise log-likelihood of model", id_model, "/ chain", i_chain, "-----\n\n\n")
  
  result <- tryCatch(
    {
      # Use withCallingHandlers to capture messages/warnings
      withCallingHandlers({
        
        cat("Starting computing\n\n")
      
        ## Set the model environment
        environment(get_summary_pointwise_model) <- models_setup[[id_model]]
        
        ## Get the model residuals and pointwise log-likelihood matrices
        out <- get_summary_pointwise_model(model_output$outputs)
        
        ## Add model IDs to output
        out$id_model <- id_model
        out$i_chain <- i_chain
        
        cat("\n\nFinished computing\n\n")
        
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
  
  
  # Return the list of models' summary matrices ----
  return(result) 

  
}