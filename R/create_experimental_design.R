create_experimental_design <- function() {
  
  # Angiosperm and gymnosperm diameter/competition interactive effect----
  # log(LAD) ~ DBH*BATOT + (1 | site/origin)
  expand.grid(
    
    filter_sensors = TRUE, # filter sensors based on preliminary analysis (converged = TRUE)
    
    site_rd_effect = TRUE, # hierarchical site into origin random effect
    origin_rd_effect = TRUE,
    
    batot_site_effect = TRUE, # site effect of total basal area
    dbh_effect = TRUE, # individual effect of diameter
    
    group_variable = "phylogeny", # can also be "species" for example
    intercept_per_group = TRUE,
    dbh_per_group = TRUE,
    
    # OPTION 1: group pooling: covariance effect between intercept and dbh slope
    # does not work with dbhXbatot interaction !!
    intercept_group_pooling = FALSE,
    dbh_group_pooling = FALSE,
    consider_covariance = FALSE,
    
    # OPTION 2: consider interaction between batot and dbh
    # does not work with group pooling, because it implies considering batot per group
    batot_per_group = TRUE,
    dbh_interaction_batot = TRUE
    
  ) %>% 

    dplyr::mutate(
      id_model = row_number(),
      n_iterations = 20000,
      n_burning = 0,
      n_subchains = 3,
      n_chains = 3
    ) %>%
    dplyr::relocate(id_model, n_iterations, n_burning, n_subchains, n_chains)
  
}
