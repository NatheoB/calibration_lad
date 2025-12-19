bind_output_light <- function(data_output_light_list) {
  
  # TREES INTERCEPTION
  output_light_trees <- dplyr::left_join(
    
    data_output_light_list$control$trees %>% 
      dplyr::mutate(lci = ifelse(is.na(lci), 1, lci)) %>% 
      dplyr::group_by(site, id_tree) %>% 
      dplyr::summarise(across(c(e, epot, lci), 
                              list(
                                "lower" = ~quantile(., 0.025),
                                "median" = ~median(.),
                                "upper" = ~quantile(., 0.975)
                              ))) %>% 
      dplyr::rename_with(~paste0(.x, "_control"), -c(site, id_tree)),
    
    data_output_light_list$models[[1]] %>%
      purrr::map(~.x$trees) %>% 
      dplyr::bind_rows(.id = "id_set") %>% 
      dplyr::mutate(lci = ifelse(is.na(lci), 1, lci)) %>% 
      dplyr::group_by(site, id_tree) %>% 
      dplyr::summarise(across(c(e, epot, lci), 
                              list(
                                "lower" = ~quantile(., 0.025),
                                "median" = ~median(.),
                                "upper" = ~quantile(., 0.975)
                              ))) %>% 
      dplyr::rename_with(~paste0(.x, "_model"), -c(site, id_tree)),
    
    by = c("site", "id_tree")
  )
  
  
  # LIGHT ON THE GROUND
  output_light_stand <- dplyr::left_join(
    
    data_output_light_list$control$stand %>% 
      dplyr::group_by(site) %>% 
      dplyr::summarise(across(c(pacl_min, pacl_mean, pacl_max), 
                              list(
                                "lower" = ~quantile(., 0.025),
                                "median" = ~median(.),
                                "upper" = ~quantile(., 0.975)
                              ))) %>% 
      dplyr::rename_with(~paste0(.x, "_control"), -site),
    
    data_output_light_list$models[[1]] %>%
      purrr::map(~.x$stand) %>% 
      dplyr::bind_rows(.id = "id_set") %>% 
      dplyr::group_by(site) %>% 
      dplyr::summarise(across(c(pacl_min, pacl_mean, pacl_max), 
                              list(
                                "lower" = ~quantile(., 0.025),
                                "median" = ~median(.),
                                "upper" = ~quantile(., 0.975)
                              ))) %>% 
      dplyr::rename_with(~paste0(.x, "_model"), -site),
    
    by = "site"
  )
  
  # Return both light at the tree-level and at the stand-level
  list(
    "trees" = output_light_trees,
    "stand" = output_light_stand
  )
}