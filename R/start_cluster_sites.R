start_cluster_sites <- function(setup_env, n_workers) {
  
  # Start the workers ----
  cl <- parallel::makePSOCKcluster(n_workers)
  
  # Load the packages on the workers ----
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(SamsaRaLight)
    NULL
  })
  
  # Send the data once (they stay in the global environment of the workers) ----
  parallel::clusterExport(cl, 
                          c("data_stands", "data_rad", "sl_geometries", "use_sl_geometry",
                            "data_sites_calib", "site_names", "n_threads"), 
                          envir = setup_env)
  
  # Send the functions with the global environment of the workers as enclosure ----
  # (otherwise the whole setup environment would be sent with them)
  for (fun_name in c("run_sl_standXlad", "compute_pacl_residuals_site")) {
    fun <- get(fun_name, envir = setup_env)
    environment(fun) <- globalenv()
    parallel::clusterCall(cl, assign, fun_name, fun, envir = globalenv())
  }
  
  return(cl)
}
