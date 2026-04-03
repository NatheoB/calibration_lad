compute_rci_trees <- function(sl_stand,
                              rad_site) {
  
  tmp_out <- SamsaRaLight::run_sl(
    sl_stand, 
    rad_site, 
    detailed_output = TRUE, 
    parallel_mode = TRUE
  )
  
  vect_n_tot <- rowSums(tmp_out$output$interceptions$cells$number$total)
  df_n_tot <- data.frame(
    id_tree = as.integer(names(vect_n_tot)),
    n_tot = unname(vect_n_tot)
  )
  
  vect_n_unobstructed <- rowSums(tmp_out$output$interceptions$cells$number_unobstructed$total)
  df_n_unobstructed <- data.frame(
    id_tree = as.integer(names(vect_n_unobstructed)),
    n_unobstructed = unname(vect_n_unobstructed)
  )
  
  sl_stand$trees %>% 
    dplyr::left_join(df_n_tot, by = "id_tree") %>% 
    dplyr::left_join(df_n_unobstructed, by = "id_tree") %>% 
    dplyr::mutate(rci = 1 - n_unobstructed / n_tot) %>% 
    dplyr::select(-n_tot, -n_unobstructed)
  
}
