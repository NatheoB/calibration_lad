get_sensors_pacl_sitespecificLAD <- function(lad_values,
                                             data_stands, 
                                             data_rad)
{

  # For each site and each LAD value
  site_names <- names(data_stands)
  
  pb <- txtProgressBar(min = 0, max = length(site_names)*length(lad_values), 
                       style = 3, width = 50, char = "=")
  i_pb <- 0
  
  out_pacl_sites <- setNames(vector("list", length(site_names)), site_names)
  for (site in site_names) {
    
    out_pacl_sites_lads <- setNames(vector("list", length(lad_values)), lad_values) 
    for (lad in lad_values) {
      
      # Copy the stand
      tmp_stand <- data_stands[[site]]
      
      # Set the site specific tree crown LAD 
      tmp_stand$trees$crown_lad <- lad
      
      # Run SamsaraLight
      tmp_out_sl <- 
        SamsaRaLight::run_sl(tmp_stand, 
                             data_rad[[site]],
                             sensors_only = TRUE,
                             use_torus = TRUE,
                             turbid_medium = TRUE,
                             detailed_output = FALSE,
                             parallel_mode = TRUE,
                             n_threads = NULL,
                             verbose = FALSE)
      
      # Sl sensors output for the plot
      out_pacl_sites_lads[[as.character(lad)]] <- tmp_out_sl$output$light$sensors %>% 
        dplyr::select(id_sensor, pacl)
      
      # Update progress bar
      i_pb <- i_pb + 1
      setTxtProgressBar(pb, i_pb)
    }
    
    out_pacl_sites[[site]] <- out_pacl_sites_lads %>% 
      dplyr::bind_rows(.id = "lad") %>% 
      dplyr::mutate(lad = as.numeric(lad))
    
  }
  close(pb)
  
  out_pacl_sites %>% 
    dplyr::bind_rows(.id = "site")
}