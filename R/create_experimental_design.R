create_experimental_design <- function() {
  
  # 1. Evaluate the effect of sensors filtering with simple model ----
  exp1 <- expand.grid(
    filter_sensors = c(FALSE, TRUE),
    species_specific = FALSE,
    site_rd_effect = FALSE,
    origin_rd_effect = FALSE,
    dbh_effect = FALSE,
    compet_effect = FALSE
  )  %>% 
    dplyr::mutate(
      id_exp = 1,
      n_chains = 1,
      n_iterations = 5000,
      id_mod = row_number(),
    ) %>%
    dplyr::relocate(id_exp, n_chains, n_iterations, id_mod)
  
  
  # 2. Simple model (random effect and species-specific) ----
  exp2 <- expand.grid(
    filter_sensors = TRUE,
    species_specific = c(FALSE, TRUE),
    site_rd_effect = c(FALSE, TRUE),
    origin_rd_effect = c(FALSE, TRUE),
    dbh_effect = FALSE,
    compet_effect = FALSE
  )  %>% 
    dplyr::mutate(
      id_exp = 2,
      n_chains = 1,
      n_iterations = 20000,
      id_mod = row_number(),
    ) %>%
    dplyr::relocate(id_exp, n_chains, n_iterations, id_mod)
  
  
  # FINAL EXPERIMENTAL DESIGN ____
  dplyr::bind_rows(exp1, exp2) %>% 
    dplyr::mutate(
      id_simu = paste0(id_exp, "_", id_mod)
    ) %>% 
    dplyr::relocate(id_simu)
  
}
