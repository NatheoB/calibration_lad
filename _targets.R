# _targets.R file
library(targets)

library(future)
library(future.callr)
plan(callr)

# Install SmasaraLight
# install.packages("SamsaRaLight_1.0.tar.gz",
#                  repos = NULL, type = "source")

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
  
  tar_target(CHAINS, 1),
  
  tar_target(LAD_CONTROL, 0.5),
  
  
  
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
  
  ## Preliminary analysis ----
  # simple minimization of residuals

  tar_target(lads_method1, seq(0.01, 5, by = 0.01)),
  
  tar_target(output_pacl_method1, get_sensors_pacl_sitespecificLAD(lads_method1,
                                                                   data_calib,
                                                                   data_rad,
                                                                   init_db$plots)),
  
  tar_target(output_lad_method1, fit_lad_method1(output_pacl_method1,
                                                 max(lads_method1),
                                                 data_sensors_punobs,
                                                 "output/residuals_sensors")),
  
  tar_target(best_lad_method1, mean(output_lad_method1$best_lad[output_lad_method1$converged])),
  
  
  ## Bayesian calibration ----
  
  ### Get the species to calibrate ----
  tar_target(sp_calib_occ, get_occurences_species2calib(init_db)),
  tar_target(species2calib, as.character(unique(sp_calib_occ$species))),
  
  
  ### Create the experimental design ----
  tar_target(exp_design, create_experimental_design()),
  tar_target(id_models, exp_design$id_model),

  
  ### Initialise the Bayesian setups ----
  # 465/1121 sensors had been removed
  # 3 sites had been removed : Cloture11, Cloture15, Cloture2
  tar_target(models_setup, initialise_models(exp_design,
                                             init_db$sensors,
                                             init_db$plots,
                                             data_calib,
                                             data_rad,
                                             output_lad_method1,
                                             species2calib,
                                             prior_lad = 0.5)),


  ### Run the MCMC ----
  tar_target(models_output_list, calibrate_models(models_setup,
                                                  sampling_algo = "DREAMzs",
                                                  id_model = id_models,
                                                  i_chain = CHAINS,
                                                  logs_folder = "logs"),
             pattern = cross(id_models, CHAINS),
             iteration = "list"),
  

  
  # COMPARE AND EVALUATE THE MODELS ----

  ## Compute pointwise likelihoods ----
  tar_target(models_summary_pointwise_list, get_summary_pointwise_models(models_setup,
                                                                         models_output_list,
                                                                         logs_folder = "logs"),
             pattern = map(models_output_list),
             iteration = "list"),


  
  ## Compare models with LOO-CV and WAIC ----
  tar_target(models_comparison, compare_models(models_summary_pointwise_list)),


  ## Evaluate models with RMSE ----
  tar_target(models_evaluation, evaluate_models(models_summary_pointwise_list)),
  
  
  
  
  
  # SAVE OUTPUTS ----
  tar_target(models_output_fp, save_models(exp_design,
                                           models_setup,
                                           models_output_list,
                                           models_summary_pointwise_list,
                                           models_comparison,
                                           models_evaluation,
                                           "output/calib/",
                                           "out_20251218_dbhXbatot.Rdata")),
  
  
  # COMPUTE OUTPUT VARIABLES ----
  
  ## Create parameters table ----
  tar_target(output_params, get_output_params(models_output_list,
                                              n_analysis = 15000,
                                              thinning = 15)),
  
  
  ## Compute tree-level variables ---- 
  tar_target(data_output_tree, compute_output_tree(data_calib,
                                                   models_setup,
                                                   output_params,
                                                   LAD_CONTROL)),
  
  
  ## Compute stand-level variables ---- 
  tar_target(data_output_stand, compute_output_stand(data_calib,
                                                     models_setup,
                                                     output_params,
                                                     LAD_CONTROL)),
  
  
  ## Apply SamsaRalight on output stands ----
  tar_target(data_output_light_list, compute_output_light(data_calib,
                                                          data_rad,
                                                          init_db$plots,
                                                          models_setup,
                                                          output_params,
                                                          LAD_CONTROL)),

  tar_target(data_output_light, bind_output_light(data_output_light_list)),
  
  NULL
)
