get_summary_pointwise_models <- function(models_setup,
                                         models_output,
                                         n_values) {
  
  get_summary_pointwise_model <- function(output, n_values) {
    
    # Get samples of the chains of the model (bind the chains samplings with coda argument)
    samples <- BayesianTools::getSample(output, coda = F)
    samples <- samples[(nrow(samples)-n_values+1):nrow(samples),]
    
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
  
  
  # Compute the residuals of all the last samples of each model ----
  
  # Calibrate all the models ----
  ids_simu <- names(models_setup)
  n_mods <- length(ids_simu)
  
  models_matrices <- setNames(vector("list", n_mods), ids_simu)
  for (i in 1:n_mods) {
    
    # Get the simu name
    id_simu <- ids_simu[i]
    
    # Print a message to follow the calibration process
    message(paste0("Computing residuals and pointwise log-likelihood of model ", id_simu, " - ", i, "/", n_mods, "..."))
    
    # Set the model environment
    environment(get_summary_pointwise_model) <- models_setup[[id_simu]]
    
    # Get the model residuals and pointwise log-likelihood matrices
    models_matrices[[id_simu]] <- get_summary_pointwise_model(models_output[[id_simu]]$outputs,
                                                              n_values)
  }
  
  # Return the list of models' summary matrices ----
  return(models_matrices)
  
}