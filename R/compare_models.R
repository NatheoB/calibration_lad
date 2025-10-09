compare_models <- function(models_setup,
                           models_output) {
  
  # Function to get the marginal likelihood of a model given the environment is set up ----
  compute_comparison_indicators <- function(output) {
    
    # Get samples of the chains of the model (bind the chains samplings with coda argument)
    samples <- BayesianTools::getSample(output, coda = F)
    
    # Get the number of observation (i.e. total number of sesnors over all the sites)
    n_sensors <- nrow(data_sensors %>% dplyr::bind_rows())
    
    # Get the matrix of pointwise log-likelihood for all samples of parameters
    pointwise_ll <- matrix(NA, nrow = nrow(samples), ncol = n_sensors)
    
    pb <- txtProgressBar(min = 0, max = nrow(samples),
                         style = 3, width = 50, char = "=")
    
    for (i in 1:nrow(samples)) {
      
      pointwise_ll[i,] <- compute_log_posterior(samples[i,], pointwise = TRUE, print.pb = FALSE)

      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Compute LOOC-CV and WAIC from the matrix of log-likelihood
    waic_result <- tryCatch({
      loo::waic(pointwise_ll)
    }, error = function(e) {
      message(e)
      NULL
    })
    
    loo_result <- tryCatch({
      loo::loo(pointwise_ll)
    }, error = function(e) {
      message(e)
      NULL
    })
    
    
    return(list(
      "WAIC" = waic_result,
      "LOO-CV" = loo_result
    ))
  }
  
  
  # Compute the LOO-CV and the WAIC of each model ----
  n_mods <- length(models_setup)
  models_comp <- vector("list", n_mods)
  for (i in 1:n_mods) {
    
    # Print a message to follow the calibration process
    message(paste0("Computing comparison indicators of model ", i, "/", n_mods, "..."))
    
    # Set the model environment
    environment(compute_comparison_indicators) <- models_setup[[i]]
    
    # Get the model marginal likelihood
    models_comp[[i]] <- compute_comparison_indicators(models_output[[i]])
  }
  
  
  # Return the list of models' output and settings ----
  return(models_comp)
}