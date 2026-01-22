compute_output_light <- function(data_stands,
                                 data_rad,
                                 models_setup,
                                 output_params,
                                 lad_control) {
  
  # Compute light for the control LAD
  message("Computing light with control LAD...")
  
  output_control <- data_stands %>% 
    purrr::map2(names(data_stands), function(sl_stand, site_name) {
      
      # Get the trees and set all trees to the control LAD
      sl_stand$trees <- sl_stand$trees %>% 
        dplyr::mutate(crown_lad = lad_control)
      
      # Run SamsaraLight
      tmp_sl_out <- SamsaRaLight::run_sl(sl_stand, 
                                         data_rad[[site_name]],
                                         sensors_only = FALSE,
                                         turbid_medium = TRUE,
                                         use_torus = TRUE,
                                         detailed_output = FALSE,
                                         parallel_mode = TRUE,
                                         n_threads = NULL,
                                         verbose = FALSE)
      
      list(
        "trees" = tmp_sl_out$output$light$trees %>% 
          dplyr::filter(id_tree %in% sl_stand$trees$id_tree[!sl_stand$trees$added_to_fill]) %>% 
          dplyr::select(id_tree, e, epot, lci),
        
        "stand" = tmp_sl_out$output$light$cells %>%
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
      
      message("Computing light for model ", id_model, "...")
      
      ids_set <- unique(out_model$id_set)
      # ids_set <- sample(ids_set, 200)
      
      output_light_sets <- setNames(vector("list", length(ids_set)), ids_set)
      for (i in seq_along(ids_set)) {
        
        print(paste("set", i, "/", length(ids_set)))
        
        tmp_output_light <- data_stands %>% 
          purrr::map2(names(data_stands), function(sl_stand, site_name) {
            
            # Get param set 
            tmp_params <- out_model[out_model$id_set == ids_set[i],] %>% 
              tidyr::pivot_longer(contains(c("intercept.", "dbh.", "batot.", "dbhXbatot.")),
                                  names_pattern = "(.*)\\.(.*)",
                                  names_to = c(".value", "phylogeny")) %>% 
              dplyr::rename_with(~paste0("p_", .), all_of(c("intercept", "dbh", "batot", "dbhXbatot")))
            
            
            # Get the trees and set all trees to the control LAD
            sl_stand$trees <- sl_stand$trees %>% 
              
              dplyr::left_join(tmp_params, by = "phylogeny") %>% 
              
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
                crown_lad = exp(p_intercept + 
                                  dbh_std * p_dbh + 
                                  batot_std * p_batot +
                                  dbh_std * batot_std * p_dbhXbatot + 
                                  0.5*variance_sum)
              )
            
            # Run SamsaraLight
            tmp_sl_out <- SamsaRaLight::run_sl(sl_stand, 
                                               data_rad[[site_name]],
                                               sensors_only = FALSE,
                                               turbid_medium = TRUE,
                                               use_torus = TRUE,
                                               detailed_output = FALSE,
                                               parallel_mode = TRUE,
                                               n_threads = NULL,
                                               verbose = FALSE)
            
            list(
              "trees" = tmp_sl_out$output$light$trees %>% 
                dplyr::filter(id_tree %in% sl_stand$trees$id_tree[!sl_stand$trees$added_to_fill]) %>% 
                dplyr::select(id_tree, e, epot, lci),
              
              "stand" = tmp_sl_out$output$light$cells %>%
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