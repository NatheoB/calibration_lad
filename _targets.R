# _targets.R file
library(targets)

# library(future)
# library(future.callr)
# plan(callr)

# Source functions in R folder
lapply(grep("R$", list.files("R", recursive = TRUE), value = TRUE), function(x) source(file.path("R", x)))

# Set options (i.e. clustermq.scheduler for multiprocess computing)
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")

tar_option_set(packages = c("dplyr", "tidyr", "data.table", 
                            "vroom", "purrr", 
                            "ggplot2", "ggforce", "sf",
                            "SamsaRaLight",
                            "BayesianTools", "extraDistr", "loo",
                            "rlang"))

# List of targets
list(
  
  # PARAMETERS ----
  
  ## Global parameters ----
  tar_target(SEED, 5030),
  
  ## SL parameters ----

  
  ## Experimental parameters ----

  
  ## Calibration parameters ----
  tar_target(n_chains, 3),
  tar_target(n_iterations, 50000),
  tar_target(n_burning, 0),

  
  
  
  # PREPARE THE CALIBRATION DATASET ----
  
  ## Load and clean the initial database ----
  tar_target(init_db_fp, "data/dataBase.RData", format = "file"),
  tar_target(plot_infos_fp, "data/plot_infos.csv", format = "file"),
  
  tar_target(init_db, load_database_from_rdata(init_db_fp,
                                               plot_infos_fp)),
  
  
  ## Set the species to calibrate ----
  tar_target(sp_calib_occ, get_occurences_species2calib(init_db)),
  tar_target(species2calib, as.character(unique(sp_calib_occ$species))),
  
  
  ## Create calibration plots from tree inventories ----
  tar_target(data_calib, create_calibration_stands(init_db,
                                                   "output/initial_sites",
                                                   SEED)),
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(init_db$plots)),
  
  

  # CALIBRATE THE LAD ----

  ## Initialise the Bayesian setups ----
  tar_target(models_setup, initialise_models(init_db$sensors,
                                             data_calib,
                                             data_rad,
                                             init_db$plots,
                                             species2calib)),
  
  ## Run the MCMC ----
  tar_target(models_output, calibrate_models(models_setup$setups,
                                             n_chains,
                                             n_iterations,
                                             n_burning,
                                             sampling_algo = "DREAMzs")),
  
  
  ## Compare models with LOO-CV and WAIC ----
  # tar_target(models_comparison, compare_models(models_setup$setups,
  #                                              models_output)),
  
  
  ## Evaluate models
  
  NULL
  )
