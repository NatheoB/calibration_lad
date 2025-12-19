create_experimental_design <- function() {
  
  # Species-specific intercept, dbh and compet effect without covariance ----
  # log(LAD) ~ DBH*BATOT + (1 | site/origin)
  expand.grid(
    
    filter_sensors = TRUE,
    
    site_rd_effect = TRUE,
    origin_rd_effect = TRUE,
    batot_in_site_effect = TRUE,
    
    intercept_per_sp = FALSE,
    intercept_sp_pooling = FALSE,
    
    dbh_effect = TRUE,
    dbh_interaction_batot = TRUE,
    dbh_per_sp = FALSE,
    dbh_sp_pooling = FALSE,
    
    consider_covariance = FALSE,
    
    compet_effect = FALSE
    
  ) %>% 

    dplyr::mutate(
      id_model = row_number(),
      n_iterations = 100000,
      n_analysis = 5000
    ) %>%
    dplyr::relocate(id_model, n_iterations, n_analysis)
  
}
