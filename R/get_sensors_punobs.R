get_sensors_punobs <- function(data_stands, 
                               data_rad)
{
  
  # For each site and each LAD value
  site_names <- names(data_stands)
  
  out_punobs_sites <- setNames(vector("list", length(site_names)), site_names)
  for (site in site_names) {
    
    # Run SamsaraLight
    tmp_out_sl <- 
      SamsaRaLight::run_sl(data_stands[[site]], 
                           data_rad[[site]],
                           sensors_only = TRUE,
                           use_torus = TRUE,
                           turbid_medium = TRUE,
                           detailed_output = FALSE,
                           parallel_mode = TRUE,
                           n_threads = NULL,
                           verbose = FALSE)
    
    # Sl sensors punobs for the plot
    # i.e. percentage of pacl from unobstructed rays
    out_punobs_sites[[site]] <- tmp_out_sl$output$light$sensors %>% 
      dplyr::select(id_sensor, punobs)
  }
  
  out_punobs_sites %>% 
    dplyr::bind_rows(.id = "site")
}