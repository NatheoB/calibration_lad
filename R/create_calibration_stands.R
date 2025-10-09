create_calibration_stands <- function(init_db, 
                                      output_plots_fp,
                                      seed) {
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Set the cell size, which will not impact the study as we compute light only on sensors
  cell_size <- 5
  
  # Output list
  out_stands <- setNames(
    vector("list", length = nrow(init_db$plots)),
    init_db$plots$name
  )
  
  
  # Initializes the progress bar
  i <- 0
  pb <- txtProgressBar(min = 0, max = nrow(init_db$plots), 
                       style = 3, width = 50, char = "=")   
  
  # For each site 
  for (site in names(out_stands)) {
    
    out_stands[[site]] <- tryCatch({
      
      # Get trees inventory
      tmp_data_trees <- init_db$trees[[site]] %>% 
        dplyr::left_join(init_db$species[[site]], 
                         by = c("SpCode" = "SpCode_SamsaraLL")) %>% 
        dplyr::mutate(crown_type = "8E") %>% 
        dplyr::select(id_tree = Id, 
                      species_code = SpCode,
                      species_calib = sp_calib,
                      x = X, y = Y, dbh_cm = Dbh,
                      crown_type, h_m = H, hbase_m = CBH, hmax_m = CMRH,
                      rn_m = RN, rs_m = RS, re_m = RE, rw_m = RW,
                      crown_openess = CrownOpenness, crown_lad = LAD)
      
      
      # If the plot extent dataframe is empty, set it to the default one
      tmp_plot_extent <- init_db$plot_extents[[site]]
      
      if (!is.null(tmp_plot_extent)) {
        tmp_plot_extent <- tmp_plot_extent %>% 
          dplyr::select(x = X, y = Y)
      }
      
      use_rect_zone <- ifelse(is.null(tmp_plot_extent), TRUE, FALSE)
      
      
      ### ONLY FOR THE SAKE OF GRAPHICAL REPRESENTATION ####
      
      # Create a square plot with the inventory in the center
      tmp_stand_unfilled <- SamsaRaLight::create_rect_stand(tmp_data_trees, 
                                                            cell_size, 
                                                            tmp_plot_extent,
                                                            use_rect_zone,
                                                            fill_around = FALSE)
      
      
      # Add sensors and shift coordinates inside the rect stand
      tmp_stand_unfilled$sensors <- init_db$sensors[[site]] %>% 
        dplyr::select(id_sensor = id, x, y, h_m = z) %>% 
        dplyr::mutate(
          x = x + tmp_stand_unfilled$info$shift_x,
          y = y + tmp_stand_unfilled$info$shift_y
        )
      
      # Plot the stand and save into the output folder as png image
      plot_stand(tmp_stand_unfilled, init_db$species[[site]], 
                 output_plots_fp, 
                 paste0(site, "_inv"))
      
      
      
      ### SAMSARALIGHT VIRTUAL PLOTS ###
      
      # Output list containing all replicates of a site
      
      # Create a square plot with the inventory in the center
      # And fill with trees around the core polygon
      tmp_stand_filled <- SamsaRaLight::create_rect_stand(tmp_data_trees, 
                                                          cell_size, 
                                                          tmp_plot_extent,
                                                          use_rect_zone,
                                                          fill_around = TRUE)
      
      # Add sensors and shift coordinates inside the rect stand
      tmp_stand_filled$sensors <- init_db$sensors[[site]] %>% 
        dplyr::select(id_sensor = id, x, y, h_m = z) %>% 
        dplyr::mutate(
          x = x + tmp_stand_filled$info$shift_x,
          y = y + tmp_stand_filled$info$shift_y
        )
      
      # Plot the stand and save into the output folder as png image
      plot_stand(tmp_stand_filled, init_db$species[[site]], 
                 output_plots_fp, 
                 paste0(site, "_virtualstand"))
      
      
      # Compute competition variables
      tmp_stand_filled$trees <- compute_competition_trees(tmp_stand_filled$trees,
                                                          tmp_stand_filled$info$new_area_ha)
      
      
      
      # Add the virtual plot to the list
      tmp_stand_filled
      
    }, error = function(e) {
      # If ther is an error, stop the pipeline
      message(paste("ERROR in site", site))
      stop(e)
    })
    
    # Update the progress bar
    i <- i+1
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  out_stands
}
