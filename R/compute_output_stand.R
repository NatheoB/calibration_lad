compute_output_stand <- function(data_stands, 
                                 models_setup,
                                 output_params,
                                 lad_m2m3_control) {
  
  # Site names
  site_names <- names(data_stands)
  
  # Parameters name of all models
  par_names <- unique(unlist(
    output_params %>% 
      purrr::map(~names(.x))
  ))
  
  pars <- setNames(rep(0, length(par_names)), par_names)
  
  
  # Compute tree volumes 
  data_volume <- data_stands %>%
    purrr::map(~.x$trees %>% dplyr::filter(!added_to_fill)) %>% 
    dplyr::bind_rows(.id = "site_name") %>% 
    dplyr::mutate(
      # Volumes of the upper fourths
      c_top = h_m - hmax_m,
      v_upper_n = (4/3) * pi * rn_m^2 * c_top / 8,
      v_upper_s = (4/3) * pi * rs_m^2 * c_top / 8,
      v_upper_w = (4/3) * pi * rw_m^2 * c_top / 8,
      v_upper_e = (4/3) * pi * re_m^2 * c_top / 8,
      v_upper = v_upper_n + v_upper_s + v_upper_w + v_upper_e,
      
      # Volumes of the lower fourths
      c_bot = hmax_m - hbase_m,
      v_lower_n = (4/3) * pi * rn_m^2 * c_bot / 8,
      v_lower_s = (4/3) * pi * rs_m^2 * c_bot / 8,
      v_lower_w = (4/3) * pi * rw_m^2 * c_bot / 8,
      v_lower_e = (4/3) * pi * re_m^2 * c_bot / 8,
      v_lower = v_lower_n + v_lower_s + v_lower_w + v_lower_e,
      
      # Volume of the assymetric crown
      volume_m3 = v_upper + v_lower
    ) %>% 
    dplyr::select(site_name, id_tree, phylogeny, species, batot_m2ha, dbh_cm, volume_m3)
  
  
  # For each model
  output_params %>% 
    purrr::map2(names(output_params), function(out_model, id_model) {
      
      out_sites <- setNames(vector("list", length(site_names)), site_names)
      
      for (site in site_names) {
        
        message("Creating output stands for model ", id_model, " in ", site, "...")
        
        
        # Add parameters of other models and set it to 0
        out_model <- tibble::add_column(out_model, !!!pars[setdiff(names(pars), names(out_model))])
        
        # Compute tree leaf area density and tree leaf area
        out_sites[[site]] <- tidyr::crossing(
          data_volume %>% dplyr::filter(site_name == site), 
          out_model
        ) %>% 
          
          dplyr::mutate(
            
            # Exponentiate sigmas
            sigma_mod = exp(sigma_log),
            sigma_site = exp(sigma_site_log),
            sigma_origin = exp(sigma_origin_log),
            
            # Compute variance (do not forget that we fit a log_normal LAD, thus consider variance term when back transforming)
            variance_sum = sigma_mod^2 + sigma_site^2 + sigma_origin^2,
            
            # Restandardize predictors
            dbh_std = (dbh_cm - models_setup[[id_model]]$dbh_mean) / models_setup[[id_model]]$dbh_sd,
            batot_std = (batot_m2ha - models_setup[[id_model]]$batot_mean) / models_setup[[id_model]]$batot_sd,
            
          ) %>% 
          
          tidyr::pivot_longer(contains(c("intercept.", "dbh.", "batot.", "dbhxbatot.")),
                              names_pattern = "(.*)\\.(.*)",
                              names_to = c(".value", "phylogeny_params")) %>% 
          
          dplyr::filter(phylogeny_params != phylogeny) %>% 
          
          dplyr::mutate(
            lad_m2m3_model = exp(intercept + 
                                   dbh_std * dbh + 
                                   batot_std * batot +
                                   dbh_std * batot_std * dbhXbatot + 
                                   0.5*variance_sum),
            
            la_m2_model = volume_m3 * lad_m2m3_model,
            la_m2_control = volume_m3 * lad_m2m3_control
          ) %>% 
          
          # COMPUTE STAND-LEVEL VARIABLES (LATOT AND LAI)
          dplyr::group_by(site_name, batot_m2ha, rep, id_params) %>% 
          dplyr::summarise(
            latot_m2_model = sum(la_m2_model),
            latot_m2_control = sum(la_m2_control)
          ) %>% 
          dplyr::ungroup() %>% 
          
          dplyr::mutate(
            area_m2 = data_stands[[site]]$transform$core_area_ha * 10000,
            lai_m2m2_model = latot_m2_model / area_m2,
            lai_m2m2_control = latot_m2_control / area_m2,
            diff_lai_m2m2 = lai_m2m2_control - lai_m2m2_model
          ) %>% 
          
          # Compute summaries for LAI and LA
          dplyr::group_by(site_name, batot_m2ha) %>% 
          dplyr::summarise(
            
            lai_model_lower = quantile(lai_m2m2_model, 0.025),
            lai_model_median = median(lai_m2m2_model),
            lai_model_upper = quantile(lai_m2m2_model, 0.975),
            
            lai_control_lower = quantile(lai_m2m2_control, 0.025),
            lai_control_median = median(lai_m2m2_control),
            lai_control_upper = quantile(lai_m2m2_control, 0.975),
            
            diff_lai_lower = quantile(diff_lai_m2m2, 0.025),
            diff_lai_median = median(diff_lai_m2m2),
            diff_lai_upper = quantile(diff_lai_m2m2, 0.975)
            
          ) %>% 
          dplyr::ungroup()
        
      }
      
      out_sites %>% dplyr:: bind_rows()
      
    }) %>% 
    
    # Bind chains
    dplyr::bind_rows()
  
  
}
