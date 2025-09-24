compute_residuals <- function(data_sensors, out_sl) {
  
  # Init output list
  out_residuals <- vector("list", length = length(out_sl))
  
  for (i in 1:length(out_sl)) {
    
    # Get estimated pacl from virtual sensors
    tmp_field_sensors <- out_sl[[i]]$sensor %>%
      dplyr::select(id_sensor, pacl_sl = pacl_slope)
    
    # Get measured pacl from field sensor
    tmp_sl_sensors <- data_sensors[[out_sl[[i]]$site]] %>%
      dplyr::select(id_sensor = id, pacl_field = PACLtotal)
    
    # Compute mean residuals between all sensors
    tmp_res_pacl <- tmp_field_sensors %>% 
      dplyr::left_join(tmp_sl_sensors, by = "id_sensor") %>% 
      dplyr::mutate(res_pacl = pacl_sl - pacl_field) %>% 
      dplyr::summarise(res_pacl = mean(res_pacl)) %>% 
      dplyr::pull(res_pacl)
         
    # Output the mean residual
    out_residuals[[i]] <- data.frame(
      site = out_sl[[i]]$site,
      replicate = out_sl[[i]]$replicate,
      lad = out_sl[[i]]$lad,
      res_pacl = tmp_res_pacl
    )
  }
  
  dplyr::bind_rows(out_residuals)
}