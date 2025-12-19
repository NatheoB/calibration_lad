compute_output_tree <- function(data_calib, 
                                models_setup,
                                output_params,
                                lad_m2m3_control) {
  
  # Parameters name of all models
  par_names <- unique(unlist(
    output_params %>% 
      purrr::map(~names(.x))
  ))
  
  pars <- setNames(rep(0, length(par_names)), par_names)
  
  
  # Compute tree volumes 
  data_volume <- data_calib %>%
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
    dplyr::select(site_name, id_tree, species_calib, batot_m2ha, dbh_cm, volume_m3)
  
  
  # For each model
  output_params %>% 
    purrr::map2(names(output_params), function(out_model, id_model) {
      
      # Add parameters of other models and set it to 0
      out_model <- tibble::add_column(out_model, !!!pars[setdiff(names(pars), names(out_model))])
      
      # Compute tree leaf area density and tree leaf area
      tidyr::crossing(
        data_volume, 
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
        
        dplyr::mutate(
          lad_m2m3_model = exp(intercept + 
                           dbh_std * dbh + 
                           batot_std * batot +
                           dbh_std * batot_std * dbhXbatot + 
                           0.5*variance_sum),
          
          la_m2_model = volume_m3 * lad_m2m3_model,
          la_m2_control = volume_m3 * lad_m2m3_control,
          
          diff_lad_m2m3 = lad_m2m3_control - lad_m2m3_model,
          diff_la_m2 = la_m2_control - la_m2_model
        )  %>% 
        
        dplyr::group_by(chain, site_name, id_tree, species_calib, 
                        batot_m2ha, dbh_cm, volume_m3) %>% 
        dplyr::summarise(
          
          # LAD summaries
          lad_model_lower = quantile(lad_m2m3_model, 0.025),
          lad_model_median = median(lad_m2m3_model),
          lad_model_upper = quantile(lad_m2m3_model, 0.975),
          
          lad_control_lower = quantile(lad_m2m3_control, 0.025),
          lad_control_median = median(lad_m2m3_control),
          lad_control_upper = quantile(lad_m2m3_control, 0.975),
          
          # Leaf area summaries
          la_model_lower = quantile(la_m2_model, 0.025),
          la_model_median = median(la_m2_model),
          la_model_upper = quantile(la_m2_model, 0.975),
          
          la_control_lower = quantile(la_m2_control, 0.025),
          la_control_median = median(la_m2_control),
          la_control_upper = quantile(la_m2_control, 0.975),
          
          # Diff between LAD/LA control and predicted
          diff_lad_lower = quantile(diff_lad_m2m3, 0.025),
          diff_lad_median = median(diff_lad_m2m3),
          diff_lad_upper = quantile(diff_lad_m2m3, 0.975),
          
          diff_la_lower = quantile(diff_la_m2, 0.025),
          diff_la_median = median(diff_la_m2),
          diff_la_upper = quantile(diff_la_m2, 0.975)
          
        ) %>% 
        
        dplyr::mutate(
          order = case_match(species_calib,
                             c("Abies_alba", "Larix_decidua",
                               "Picea_abies", "Pinus_sylvestris",
                               "Pseudotsuga_menziesii") ~ "gymnosperm",
                             c("Carpinus_betulus", "Fagus_sylvatica",
                               "Quercus_sp", "other") ~ "angiosperm")
        ) %>% 
        dplyr::ungroup()
    }) %>% 
    
    # Bind chains
    dplyr::bind_rows()
  
  
}