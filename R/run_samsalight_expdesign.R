run_samsalight_expdesign <- function(invs_for_samsalight,
                                     data_sensors, data_rad) {
  
  inv_ids <- names(invs_for_samsalight)
  out_list <- setNames(vector("list", length(inv_ids)), inv_ids)
  
  for (id in inv_ids) {
    
    tmp_out_samsalight <- 
      sl_run(invs_for_samsalight[[id]], 
             data_rad,
             sensors = data_sensors, sensors_only = TRUE,
             latitude = 46, slope = 6, 
             aspect = 144, north_to_x_cw = 54,
             start_day = 121, end_day = 273,
             cell_size = 10, n_cells_x = 10, n_cells_y = 10,
             turbid_medium = TRUE,
             trunk_interception = FALSE,
             soc = TRUE,
             height_anglemin = 15,
             direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
             direct_anglestep = 5,
             diffuse_anglestep = 15)
    
    out_list[[id]] <- tmp_out_samsalight$sensors
    
  }
  
  out_list
}