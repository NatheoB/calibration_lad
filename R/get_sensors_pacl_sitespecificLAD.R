get_sensors_pacl_sitespecificLAD <- function(lad_values,
                                             data_calib, 
                                             data_rad, 
                                             data_plots)
{

  # For each site and each LAD value
  site_names <- data_plots$name
  
  pb <- txtProgressBar(min = 0, max = length(site_names)*length(lad_values), 
                       style = 3, width = 50, char = "=")
  i_pb <- 0
  
  out_pacl_sites <- setNames(vector("list", length(site_names)), site_names)
  for (site in site_names) {
    
    out_pacl_sites_lads <- setNames(vector("list", length(lad_values)), lad_values) 
    for (lad in lad_values) {
      
      # Set the site specific tree crown LAD 
      tmp_trees <- data_calib[[site]]$trees %>% 
        dplyr::mutate(crown_lad = lad)
      
      # Get the plot info
      tmp_plot <- data_plots %>% 
        dplyr::filter(name == site)
      
      # Run SamsaraLight
      tmp_out_samsalight <- 
        sl_run(tmp_trees, 
               data_rad[[site]],
               sensors = data_calib[[site]]$sensors, 
               sensors_only = TRUE,
               latitude = tmp_plot$latitude, 
               slope = tmp_plot$slope, 
               aspect = tmp_plot$aspect, 
               north_to_x_cw = tmp_plot$northToX,
               start_day = 121, 
               end_day = 273,
               cell_size = data_calib[[site]]$info$cell_size, 
               n_cells_x = data_calib[[site]]$info$n_cells_x, 
               n_cells_y = data_calib[[site]]$info$n_cells_y,
               turbid_medium = TRUE,
               trunk_interception = TRUE,
               soc = TRUE,
               height_anglemin = 15,
               direct_startoffset = 0,
               direct_anglestep = 5,
               diffuse_anglestep = 15,
               detailed_output = TRUE)
      
      # Sl sensors output for the plot
      out_pacl_sites_lads[[as.character(lad)]] <- tmp_out_samsalight$output$sensors %>% 
        dplyr::select(id_sensor,
                      pacl = pacl_horizontal,
                      pacl_direct = pacl_horizontal_direct,
                      pacl_diffuse = pacl_horizontal_diffuse)
      
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