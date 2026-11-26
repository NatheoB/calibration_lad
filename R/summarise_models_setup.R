summarise_models_setup <- function(models_setup) {
  
  # Light summary of each model setup, to interpret the parameters in the analyses
  # (the setup environments carry the stands and geometries, ~190 MB each, and are not saved)
  purrr::map(models_setup, function(setup) {
    c(
      list(
        "mod_design" = setup$mod_design,
        "par_names" = setup$bayesianSetup$names,
        "site_names" = setup$site_names,
        "origin_names" = setup$origin_names,
        "species2calib" = setup$species2calib,
        "groups2calib" = setup$groups2calib
      ),
      # Standardisation constants (compet_* only exist in models with a competition effect)
      mget(c("shadetol_mean", "shadetol_sd", "dbh_mean", "dbh_sd", "compet_mean", "compet_sd"),
           envir = setup, ifnotfound = list(NA_real_))
    )
  })
}
