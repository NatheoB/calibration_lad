create_experimental_design <- function(n_iterations_per_rep) {
  
  # (1) CONTROL: 
  #   log(LAD) ~ 1 + (1 | species) + (1 | origin/site)
  #
  # (2) AUTOECOLOGY: 
  #   log(LAD) ~ broad/needle-leaved * low/high-shadetol + (1 | species) + (1 | origin/site)
  #
  # (3) DBH EFFECT: 
  #   log(LAD) ~ DBH * broad/needle-leaved * low/high-shadetol + (1 | species) + (1 | origin/site)
  #
  # (4) INTERACTION COMPET: 
  #   log(LAD) ~ DBH * BAH * broad/needle-leaved * low/high-shadetol + (1 | species) + (1 | origin/site)
  #
  # (5) LIGHT COMPET: 
  #   log(LAD) ~ DBH * RCI * broad/needle-leaved * low/high-shadetol + (1 | species) + (1 | origin/site)
  #
  
  exp_design = data.frame(
    
    # Filter sensors based on preliminary analysis (converged = TRUE)
    filter_sensors        = TRUE, 
    
    # Nested site random intercept around origin-specific mean
    site_rd_effect        = TRUE,
    origin_rd_effect      = TRUE,
    
    # Species random effect
    species_rd_effect     = TRUE,

    # Slope effects : interactive effect between dbh and competition
    dbh_effect            = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    compet_effect         = c(FALSE, FALSE, FALSE, TRUE, TRUE),
    dbh_interaction_compet = c(FALSE, FALSE, FALSE, TRUE, TRUE),
    compet_var            = c(NA, NA, NA, "balh_m2ha", "rci"),
    
    # Effect of group: broad/needle-leaved and low/high-shadetol
    grouping_var          = c(NA, "species_group", "species_group", "species_group", "species_group"),
    dbh_per_group         = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    compet_per_group      = c(FALSE, FALSE, FALSE, TRUE, TRUE),
    interaction_per_group = c(FALSE, FALSE, FALSE, TRUE, TRUE)
    
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
