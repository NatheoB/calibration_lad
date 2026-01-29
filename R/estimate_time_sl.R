estimate_time_sl <- function(data_stands, 
                             data_rad, 
                             test_n_threads) {
  
  # For each site and each LAD value
  site_names <- names(data_stands)
  
  pb <- txtProgressBar(min = 0, max = length(site_names)*length(test_n_threads), 
                       style = 3, width = 50, char = "=")
  i_pb <- 0
  
  out_time_sites <- setNames(vector("list", length(site_names)), site_names)
  for (site in site_names) {
    
    out_time_sites_nthreads <- setNames(vector("list", length(test_n_threads)), as.character(test_n_threads)) 
    for (n in test_n_threads) {
      
      # Run SamsaraLight
      out_time <- microbenchmark::microbenchmark(
        
        SamsaRaLight::run_sl(data_stands[[site]], 
                             data_rad[[site]],
                             sensors_only = TRUE,
                             use_torus = TRUE,
                             turbid_medium = TRUE,
                             detailed_output = FALSE,
                             parallel_mode = TRUE,
                             n_threads = n,
                             verbose = FALSE)
        
      ) 
      
      # Save time
      s <- summary(out_time)
      
      out_time_sites_nthreads[[as.character(n)]] <- 
        data.frame(
          min = s$min,
          lq = s$lq,
          mean = s$mean,
          median = s$median,
          uq = s$uq,
          max = s$max,
          neval = s$neval
        ) 
      
      # Update progress bar
      i_pb <- i_pb + 1
      setTxtProgressBar(pb, i_pb)
    }
    
    out_time_sites[[site]] <- out_time_sites_nthreads %>% 
      dplyr::bind_rows(.id = "n_threads")
    
  }
  close(pb)
  
  out_time_sites %>% 
    dplyr::bind_rows(.id = "site")
  
}