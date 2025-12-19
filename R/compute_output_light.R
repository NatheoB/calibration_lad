compute_output_light <- function(data_calib,
                                 data_rad,
                                 data_plots,
                                 models_setup,
                                 output_params,
                                 lad_control) {
  
  # Compute light for the control LAD
  message("Computing light with control LAD...")
  
  output_control <- data_calib %>% 
    purrr::map2(names(data_calib), function(data_site, site_name) {
      
      # Get the trees and set all trees to the control LAD
      tmp_trees <- data_site$trees %>% 
        dplyr::mutate(crown_lad = lad_control)
      
      # Get the plot info
      tmp_plot <- data_plots %>% 
        dplyr::filter(name == site_name)
      
      # Run SamsaraLight
      tmp_out_samsalight <- 
        SamsaRaLight::sl_run(tmp_trees, 
                             data_rad[[site_name]],
                             sensors = NULL, 
                             sensors_only = FALSE,
                             latitude = tmp_plot$latitude, 
                             slope = tmp_plot$slope, 
                             aspect = tmp_plot$aspect, 
                             north_to_x_cw = tmp_plot$northToX,
                             start_day = 1, 
                             end_day = 365,
                             cell_size = data_site$info$cell_size, 
                             n_cells_x = data_site$info$n_cells_x, 
                             n_cells_y = data_site$info$n_cells_y,
                             turbid_medium = TRUE,
                             trunk_interception = TRUE,
                             soc = TRUE,
                             height_anglemin = 15,
                             direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
                             direct_anglestep = 5,
                             diffuse_anglestep = 15,
                             detailed_output = FALSE)
      
      list(
        "trees" = tmp_out_samsalight$output$trees %>% 
          dplyr::filter(id_tree %in% tmp_trees$id_tree[!tmp_trees$added_to_fill]) %>% 
          dplyr::select(id_tree, e, epot, lci),
        
        "stand" = tmp_out_samsalight$output$cells %>%
          dplyr::summarise(
            pacl_min = min(pacl),
            pacl_mean = mean(pacl),
            pacl_sd = sd(pacl),
            pacl_max = max(pacl)
          )
      )
    })
  
  
  output_control <- list(
    
    "trees" = output_control %>% 
      purrr::map_dfr(
        ~ .x$trees,
        .id = "site"
      ),
    
    "stand" = output_control %>% 
      purrr::map_dfr(
        ~ .x$stand,
        .id = "site"
      )
  )
  
  
  # Parameters name of all models
  par_names <- unique(unlist(
    output_params %>% 
      purrr::map(~names(.x))
  ))
  
  pars <- setNames(rep(0, length(par_names)), par_names)
  
  
  # For each model
  output_light <- output_params %>% 
    purrr::map2(names(output_params), function(out_model, id_model) {
      
      # Add parameters of other models and set it to 0
      out_model <- tibble::add_column(out_model, !!!pars[setdiff(names(pars), names(out_model))])
      
      message(paste("Computing light for LAD model", id_model))
      
      ids_set <- unique(out_model$id_set)
      # ids_set <- sample(ids_set, 100)
      
      output_light_sets <- setNames(vector("list", length(ids_set)), ids_set)
      for (i in seq_along(ids_set)) {
        
        print(paste("set", i, "/", length(ids_set)))
        
        tmp_output_light <- data_calib %>% 
          purrr::map2(names(data_calib), function(data_site, site_name) {
            
            # Get param set 
            tmp_params <- out_model[out_model$id_set == ids_set[i],]
            
            # Get the trees and set all trees to the control LAD
            tmp_trees <- data_site$trees %>% 
              
              dplyr::mutate(
                
                # Exponentiate sigmas
                sigma_mod = exp(tmp_params$sigma_log),
                sigma_site = exp(tmp_params$sigma_site_log),
                sigma_origin = exp(tmp_params$sigma_origin_log),
                
                # Compute variance (do not forget that we fit a log_normal LAD, thus consider variance term when back transforming)
                variance_sum = sigma_mod^2 + sigma_site^2 + sigma_origin^2,
                
                # Restandardize predictors
                dbh_std = (dbh_cm - models_setup[[id_model]]$dbh_mean) / models_setup[[id_model]]$dbh_sd,
                batot_std = (batot_m2ha - models_setup[[id_model]]$batot_mean) / models_setup[[id_model]]$batot_sd,
                
              ) %>% 
              
              dplyr::mutate(
                crown_lad = exp(tmp_params$intercept + 
                                  dbh_std * tmp_params$dbh + 
                                  batot_std * tmp_params$batot +
                                  dbh_std * batot_std * tmp_params$dbhXbatot + 
                                  0.5*variance_sum)
              )
            
            
            
            
            # Get the plot info
            tmp_plot <- data_plots %>% 
              dplyr::filter(name == site_name)
            
            # Run SamsaraLight
            tmp_out_samsalight <- 
              SamsaRaLight::sl_run(tmp_trees, 
                                   data_rad[[site_name]],
                                   sensors = NULL, 
                                   sensors_only = FALSE,
                                   latitude = tmp_plot$latitude, 
                                   slope = tmp_plot$slope, 
                                   aspect = tmp_plot$aspect, 
                                   north_to_x_cw = tmp_plot$northToX,
                                   start_day = 1, 
                                   end_day = 365,
                                   cell_size = data_site$info$cell_size, 
                                   n_cells_x = data_site$info$n_cells_x, 
                                   n_cells_y = data_site$info$n_cells_y,
                                   turbid_medium = TRUE,
                                   trunk_interception = TRUE,
                                   soc = TRUE,
                                   height_anglemin = 15,
                                   direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
                                   direct_anglestep = 5,
                                   diffuse_anglestep = 15,
                                   detailed_output = FALSE)
            
            list(
              "trees" = tmp_out_samsalight$output$trees %>% 
                dplyr::filter(id_tree %in% tmp_trees$id_tree[!tmp_trees$added_to_fill]) %>% 
                dplyr::select(id_tree, e, epot, lci),
              
              "stand" = tmp_out_samsalight$output$cells %>%
                dplyr::summarise(
                  pacl_min = min(pacl),
                  pacl_mean = mean(pacl),
                  pacl_sd = sd(pacl),
                  pacl_max = max(pacl)
                )
            )
          })
        
        output_light_sets[[ids_set[i]]] <- list(
          
          "trees" = tmp_output_light %>% 
            purrr::map_dfr(
              ~ .x$trees,
              .id = "site"
            ),
          
          "stand" = tmp_output_light %>% 
            purrr::map_dfr(
              ~ .x$stand,
              .id = "site"
            )
        )
        
      }
      
      output_light_sets
    })
  
  
  # Combine all computed lights
  list(
    "control" = output_control,
    "models" = output_light
  )
}