run_samsalight_expdesign <- function(data_sl,
                                     data_rad,
                                     exp_design) {
  
  # Output list
  out_list <- vector("list", nrow(exp_design))
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nrow(exp_design), 
                       style = 3, width = 50, char = "=")   
  
  # Run for each simulation of the experimental design
  for (i in 1:nrow(exp_design)) {
    
    # Get info about the simulation
    tmp_simu_site <- as.character(exp_design$inv_name[i])
    tmp_simu_rep <- exp_design$replicate[i]
    tmp_simu_lad <- exp_design$lad[i]
    
    # Get tree dataset and set the LAD value for all trees
    tmp_trees <- data_sl[[tmp_simu_site]][[tmp_simu_rep]]$trees %>% 
      dplyr::mutate(crown_lad = tmp_simu_lad)
    
    # Run SamsaraLight
    tmp_out_samsalight <- 
      sl_run(tmp_trees, 
             data_rad[[tmp_simu_site]],
             sensors = data_sl[[tmp_simu_site]][[tmp_simu_rep]]$sensors, 
             sensors_only = TRUE,
             latitude = 46, slope = 6, 
             aspect = 144, north_to_x_cw = 54,
             start_day = 121, end_day = 273,
             cell_size = data_sl[[tmp_simu_site]][[tmp_simu_rep]]$info$cell_size, 
             n_cells_x = data_sl[[tmp_simu_site]][[tmp_simu_rep]]$info$n_cells_x, 
             n_cells_y = data_sl[[tmp_simu_site]][[tmp_simu_rep]]$info$n_cells_y,
             turbid_medium = TRUE,
             trunk_interception = FALSE,
             soc = TRUE,
             height_anglemin = 15,
             direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
             direct_anglestep = 5,
             diffuse_anglestep = 15)
    
    # Save the sl sensor output
    out_list[[i]] <- list(
      "site" = tmp_simu_site,
      "replicate" = tmp_simu_rep,
      "lad" = tmp_simu_lad,
      "sensor" = tmp_out_samsalight$sensors
    )
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  out_list
}