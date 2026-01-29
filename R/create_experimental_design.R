create_experimental_design <- function() {
  
  # Angiosperm/gymnosperm diameter/competition interactive effect----
  # log(LAD) ~ DBH*BATOT + group + (1 | species) + (1 | site/origin)
  data.frame(
    
    # Filter sensors based on preliminary analysis (converged = TRUE)
    filter_sensors        = TRUE, 
    
    # Hierarchical site into origin random effect
    site_rd_effect        = TRUE,
    origin_rd_effect      = TRUE,

    # Slope effects dbh (individual-level), batotal (site-level) and interaction
    batot_site_effect     = TRUE,
    dbh_effect            = TRUE,
    dbh_interaction_batot = TRUE,
    
    # Species-pooling (intercept_sp_pooling is the same as a species random effect)
    intercept_sp_pooling  = TRUE,
    dbh_sp_pooling        = FALSE,
    covariance_sp_pooling = FALSE,
    
    # Effect of functional group: angio/gymno (angiosperm is the reference group)
    group_effect          = TRUE,
    dbh_per_group         = FALSE,
    batot_per_group       = FALSE,
    interaction_per_group = FALSE
    
  ) %>% 

    dplyr::mutate(
      id_model = row_number(),
      n_chains = 1,
      n_iterations = 150000,
      n_burning = 0,
      nCR = 3
    ) %>%
    dplyr::relocate(id_model, n_chains, n_iterations, n_burning, nCR)
  
}
