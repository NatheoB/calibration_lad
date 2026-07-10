compute_output_tree_light <- function(data_stands,
                                      data_rad,
                                      output_tree_lad) {
  
  id_sets <- unique(output_tree_lad$out_lad$id_set)
  
  # Estimate light with SamsaraLight
  out_light_list <- vector("list", length(id_sets))
  for (i in seq_along(id_sets)) {
    
    sl_stand <- data_stands[[output_tree_lad$site_name]]
    
    sl_stand$trees <- sl_stand$trees %>% 
      dplyr::left_join(output_tree_lad$out_lad %>% 
                         dplyr::filter(id_set == id_sets[i]),
                       by = "id_tree") %>% 
      dplyr::mutate(crown_lad = lad)
    
    out_sl <- SamsaRaLight::run_sl(sl_stand, 
                                   data_rad[[output_tree_lad$site_name]],
                                   sensors_only = FALSE,
                                   detailed_output = FALSE,
                                   parallel_mode = TRUE,
                                   n_threads = NULL,
                                   verbose = FALSE)
    
    out_light_list[[i]] <- list(
      "out_cells" = out_sl$output$light$cells %>% 
        dplyr::mutate(id_set = id_sets[i]) %>% 
        dplyr::select(id_cell, id_set, e, pacl), 
      
      "out_trees" = out_sl$output$light$trees %>%  
        dplyr::mutate(id_set = id_sets[i]) %>% 
        dplyr::select(id_tree, id_set, e, epot)
    )
  }

  out_light <- list(
    "site_name" = output_tree_lad$site_name,
    "id_mod" = output_tree_lad$id_mod,
    "out_cells" = purrr::map_dfr(out_light_list, ~ .x$out_cells),
    "out_trees" = purrr::map_dfr(out_light_list, ~ .x$out_trees)
  )
}
