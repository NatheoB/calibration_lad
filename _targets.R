# _targets.R file
library(targets)

# library(future)
# library(future.callr)
# plan(callr)

# Source functions in R folder
lapply(grep("R$", list.files("R", recursive = TRUE), value = TRUE), 
       function(x) source(file.path("R", x)))


# Set options (i.e. clustermq.scheduler for multiprocess computing)
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")

tar_option_set(packages = c("dplyr", "tidyr", "data.table", 
                            "vroom", "purrr", 
                            "ggplot2", "ggforce", "sf",
                            "SamsaRaLight",
                            "BayesianTools", "extraDistr", "loo",
                            "rlang", "e1071", "mgcv"))

# List of targets
list(
  
  # PARAMETERS ----
  tar_target(SEED, 5030),
  tar_target(n_posteriors_analysis, 1000),
  
  
  
  # PREPARE CALIBRATION ----
  
  ## Load and clean the initial database ----
  tar_target(init_db_fp, "data/dataBase.RData", format = "file"),
  tar_target(plot_infos_fp, "data/plot_infos.csv", format = "file"),
  
  tar_target(init_db, load_database_from_rdata(init_db_fp,
                                               plot_infos_fp)),
  
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(init_db$plots)),
  
  
  ## Create calibration plots from tree inventories ----
  tar_target(data_calib, create_calibration_stands(init_db,
                                                   "output/initial_sites",
                                                   SEED)),

  ## Get the sensors to calibrate ----
  tar_target(thresholds_punobs, c(0.4, 0.8, 1)),
  
  tar_target(data_sensors_punobs, get_sensors_punobs(init_db,
                                                     data_calib,
                                                     data_rad,
                                                     thresholds_punobs,
                                                     "output/initial_sensors")),
  
  # CALIBRATE THE LAD ----
  
  ## METHOD 1: simple minimization of residuals ----

  tar_target(lads_method1, seq(0.001, 5, by = 0.001)),
  
  tar_target(output_pacl_method1, get_sensors_pacl_sitespecificLAD(lads_method1,
                                                                   data_calib,
                                                                   data_rad,
                                                                   init_db$plots)),
  
  tar_target(output_lad_method1, fit_lad_method1(output_pacl_method1, 
                                                 data_sensors_punobs,
                                                 "output/residuals_sensors")),
  
  
  ## METHOD 2: Bayesian calibration ----
  
  ### Get the species to calibrate ----
  tar_target(sp_calib_occ, get_occurences_species2calib(init_db)),
  tar_target(species2calib, as.character(unique(sp_calib_occ$species))),
  
  
  ### Create the experimental design ----
  tar_target(exp_design, create_experimental_design()),


  ### Initialise the Bayesian setups ----
  # 465/1121 sensors had been removed
  # 3 sites had been removed : Cloture11, Cloture15, Cloture2
  tar_target(models_setup, initialise_models(exp_design,
                                             init_db$sensors,
                                             init_db$plots,
                                             data_calib,
                                             data_rad,
                                             output_lad_method1,
                                             species2calib)),


  ### Run the MCMC ----
  tar_target(models_output, calibrate_models(models_setup,
                                             sampling_algo = "DREAMzs")),


  # COMPARE AND EVALUATE THE MODELS ----

  ## Get pointwise matrices (log-likelihood and residuals) ----
  tar_target(models_summary_pointwise, get_summary_pointwise_models(models_setup,
                                                                    models_output,
                                                                    n_posteriors_analysis)),


  ## Compare models with LOO-CV and WAIC ----
  tar_target(models_comparison, compare_models(models_summary_pointwise)),


  ## Evaluate models with RMSE ----
  tar_target(models_evaluation, evaluate_models(models_summary_pointwise)),
  
  
  NULL
)
