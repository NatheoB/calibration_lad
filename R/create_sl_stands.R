create_sl_stands <- function(init_db,
                             cell_size,
                             seed) {
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Output list
  out_stands <- setNames(
    vector("list", length = nrow(init_db$plots)),
    init_db$plots$name
  )
  

  # For each site 
  for (site in names(out_stands)) {
    
    message("---- Creating SamsaRaLight stand for ", site, "...")
    
    # Format and check trees inventory
    tmp_inv <- init_db$trees[[site]] %>% 
      dplyr::left_join(init_db$species[[site]], 
                       by = c("SpCode" = "SpCode_SamsaraLL")) %>% 
      dplyr::mutate(crown_type = "8E") %>% 
      dplyr::select(id_tree = Id, 
                    species_in_inv = Essence_Latin,
                    functional_group,
                    species = sp_calib,
                    x = X, y = Y, dbh_cm = Dbh,
                    crown_type, h_m = H, hbase_m = CBH, hmax_m = CMRH,
                    rn_m = RN, rs_m = RS, re_m = RE, rw_m = RW,
                    crown_openness = CrownOpenness, crown_lad = LAD)
    
    SamsaRaLight::check_inventory(tmp_inv)
  
    
    
    # Set the core polygon dataframe
    tmp_plot_extent <- init_db$plot_extents[[site]]
    
    if (!is.null(tmp_plot_extent)) {
      tmp_plot_extent <- tmp_plot_extent %>% 
        dplyr::select(x = X, y = Y)
    }
    
    
    # Get plot info
    tmp_plot_info <- init_db$plots[init_db$plots$name == site, ] 
    
    
    # Format and check sensors
    tmp_sensors <- init_db$sensors[[site]] %>% 
      dplyr::rename(id_sensor = id, h_m = z)
    
    SamsaRaLight::check_sensors(tmp_sensors)
    
    
    # Create a square plot with the inventory in the center
    tmp_stand <- SamsaRaLight::create_sl_stand(
      trees_inv = tmp_inv,
      cell_size = cell_size,
      latitude = tmp_plot_info$latitude,
      slope = tmp_plot_info$slope,
      aspect = tmp_plot_info$aspect,
      north2x = tmp_plot_info$northToX,
      sensors = tmp_sensors,
      core_polygon_df = tmp_plot_extent,
      aarect_zone = tmp_plot_info$aarect_zone, # for all except cloture (circular plots)
      fill_around = TRUE
    )
    
    # Compute competition variables
    tmp_stand$trees <- compute_competition_trees(tmp_stand$trees,
                                                 tmp_stand$transform$new_area_ha)
    
    
    # Add the virtual plot to the list
    out_stands[[site]] <- tmp_stand
  }
  
  out_stands
}
