compute_output_tree_lad <- function(data_stands, 
                                    models_setup,
                                    output_params,
                                    id_mod, 
                                    site_name) {
  
  # Keep only interest  model
  params <- output_params[[id_mod]]
  setup <- models_setup[[id_mod]]
  
  # Get species-specific and leaf group parameters
  spspe_params <- params %>%
    tidyr::pivot_longer(
      starts_with("p_z_sp."),
      names_to = "species",
      values_to = "z_sp",
      names_prefix = "p_z_sp."
    ) %>%
    tidyr::pivot_longer(
      starts_with(c("p_intercept", "p_dbh", "p_shadetol", "p_compet")),
      names_to = c(".value", "leaf_group"),
      names_sep = "\\."
    )
  
  # Add missing params for not full models
  param_names <- c("p_shadetol", "p_dbh", "p_compet",
                   "p_shadetolXdbh", "p_shadetolXcompet", "p_dbhXcompet",
                   "p_shadetolXdbhXcompet")
  
  missing_params <- setdiff(param_names, names(spspe_params))
  
  if (length(missing_params) > 0) {
    spspe_params <- spspe_params %>% 
      tibble::add_column(
        !!!setNames(rep(list(0), length(missing_params)), missing_params)
      )
  }
  
  # Estimate LAD of all trees for each parameter set
  out_lad <- tidyr::crossing(data_stands[[site_name]]$trees, id_set = params$id_set) %>% 
    dplyr::mutate(lad = estimate_lad(., spspe_params, setup)) %>% 
    dplyr::select(id_tree, id_set, lad)
  
  list(
    "site_name" = site_name, 
    "id_mod" = id_mod,
    "out_lad" = out_lad
  )
}


estimate_lad <- function(trees, spspe_params, setup) {
  
  trees %>% 
    
    dplyr::left_join(spspe_params, by = c("id_set", "species", "leaf_group")) %>% 

    dplyr::mutate(
      
      shadetol_scaled = scale_shadetol(shadetol, setup),
      dbh_scaled = scale_dbh(dbh_cm, setup),
      rci_scaled = scale_compet(rci, setup),
      
      sp_effect =
        z_sp * exp(p_sigma_sp_log),
      
      lad = exp(
        p_intercept +
          sp_effect + 
          shadetol_scaled * p_shadetol +
          dbh_scaled * p_dbh +
          rci_scaled * p_compet +
          shadetol_scaled * dbh_scaled * p_shadetolXdbh +
          shadetol_scaled * rci_scaled * p_shadetolXcompet +
          dbh_scaled * rci_scaled * p_dbhXcompet +
          shadetol_scaled * dbh_scaled * rci_scaled * p_shadetolXdbhXcompet
      )
    ) %>% 
    dplyr::pull(lad)
}

scale_shadetol <- function(x, setup) {
  
  if (!setup$mod_design$shadetol_effect) return(0)
  
  mean_shadetol <- setup$shadetol_mean
  sd_shadetol <- setup$shadetol_sd
  
  (x - mean_shadetol) / sd_shadetol
}

scale_dbh <- function(x, setup) {
  
  if (!setup$mod_design$dbh_effect) return(0)
  
  mean_dbh <- setup$dbh_mean
  sd_dbh <- setup$dbh_sd
  
  (x - mean_dbh) / sd_dbh
}

scale_compet <- function(x, setup) {
  
  if (!setup$mod_design$compet_effect) return(0)
  
  mean_compet <- setup$compet_mean
  sd_compet <- setup$compet_mean
  
  (x - mean_compet) / sd_compet
}
