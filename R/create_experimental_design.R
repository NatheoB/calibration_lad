create_experimental_design <- function(n_iterations_per_rep) {
  
  # (1) SHADE TOLERANCE: 
  #   log(LAD) ~ broad/needle-leaved * SHADETOL + (1 | species) + (1 | origin/site)
  #
  # (2) DBH EFFECT: 
  #   log(LAD) ~ DBH * broad/needle-leaved * SHADETOL + (1 | species) + (1 | origin/site)
  #
  # (3) COMPET EFFECT (BATOT, BALARGER, LIGHT)
  #   log(LAD) ~ COMPET * DBH * broad/needle-leaved * SHADETOL + (1 | species) + (1 | origin/site)

  
  exp_design = data.frame(
    
    # Filter sensors based on preliminary analysis (converged = TRUE)
    filter_sensors    = TRUE, 
    
    # Nested site random intercept around origin-specific mean
    site_rd_effect    = TRUE,
    origin_rd_effect  = TRUE,
    
    # Species random effect
    species_rd_effect = TRUE,

    # Effect of leaf group: broad/needle-leaved
    grouping_var      = "leaf_group",
    
    # Slope effects : interactive effect between dbh and competition that varies with species shade tolerance
    shadetol_effect   = TRUE,
    
    dbh_effect        = c(FALSE, TRUE, TRUE, TRUE, TRUE),
    
    compet_effect     = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    compet_var        = c(NA, NA, "batot_local_m2ha", "bal_local_m2ha", "rci")
    
  ) %>% 

    dplyr::mutate(
      id_model = row_number(),
      n_chains = 1,
      n_iterations = n_iterations_per_rep,
      n_burning = 0,
      nCR = 3
    ) %>%
    dplyr::relocate(id_model, n_chains, n_iterations, n_burning, nCR)
  
}
