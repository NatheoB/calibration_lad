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
  
  tar_target(LAD_CONTROL, 0.5),
  
  
  
  # PREPARE CALIBRATION ----
  
  ## Load the initial database ----
  tar_target(init_db_fp, "data/dataBase.RData", format = "file"),
  tar_target(plot_infos_fp, "data/plot_infos.csv", format = "file"),
  
  tar_target(init_db_raw, load_database_from_rdata(init_db_fp,
                                                   plot_infos_fp)),
  
  
  ## Clean the database ----
  tar_target(init_db, clean_database(init_db_raw)),
  
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(init_db$plots)),
  
  
  ## Create SamsaRaLight stands from tree inventories ----
  tar_target(data_stands, create_sl_stands(init_db,
                                           cell_size = 10,
                                           seed = SEED)),

  ## Plot the calibration data ----
  # tar_target(plots_stands_fp, plot_stands(data_stands,
  #                                         data_species = init_db$species,
  #                                         output_plots_fp = "output/initial_sites")),
  
  
  
  # PRELIMINARY ANALYSIS ----
  # simple minimization of residuals

  tar_target(lads_method1, seq(0.01, 5, by = 0.01)),
  tar_target(lad_convergence_threshold, 3),

  
  ## Estimate PACL given site LAD ----
  tar_target(output_pacl_method1, get_sensors_pacl_sitespecificLAD(lads_method1,
                                                                   data_stands,
                                                                   data_rad)),

  ## Estimate punobs per sensor ----
  tar_target(output_punobs, get_sensors_punobs(data_stands,
                                               data_rad)),
  
  ## Optimize LAD per site ----
  tar_target(output_lad_method1, fit_lad_method1(output_pacl_method1,
                                                 lad_convergence_threshold,
                                                 data_stands,
                                                 output_punobs,
                                                 "output/residuals_sensors")),


  ## Bayesian calibration ----

  ### Create the experimental design ----
  tar_target(exp_design, create_experimental_design()),
  tar_target(id_models, exp_design$id_model),


  ### Initialise the Bayesian setups ----
  tar_target(models_setup, initialise_models(exp_design,
                                             init_db$plots,
                                             data_stands,
                                             data_rad,
                                             output_lad_method1,
                                             prior_lad = LAD_CONTROL)),


  ### Run the MCMC ----
  tar_target(N_REPS, 3),
  tar_target(reps, 1:N_REPS),
  
  tar_target(models_output_list, calibrate_models(models_setup,
                                                  id_model = id_models,
                                                  i_rep = reps,
                                                  logs_folder = "logs/calib"),
             pattern = cross(id_models, reps),
             iteration = "list"),



  # COMPARE AND EVALUATE THE MODELS ----

  ## Compute pointwise likelihoods ----
  # tar_target(models_summary_pointwise_list, get_summary_pointwise_models(models_setup,
  #                                                                        models_output_list,
  #                                                                        logs_folder = "logs/pointwise"),
  #            pattern = map(models_output_list),
  #            iteration = "list"),



  ## Compare models with LOO-CV and WAIC ----
  # tar_target(models_comparison, compare_models(models_summary_pointwise_list)),


  ## Evaluate models with RMSE ----
  # tar_target(models_evaluation, evaluate_models(models_summary_pointwise_list)),





  # SAVE OUTPUTS ----
  # tar_target(models_output_fp, save_models(exp_design,
  #                                          models_setup,
  #                                          models_output_list,
  #                                          models_summary_pointwise_list,
  #                                          models_comparison,
  #                                          models_evaluation,
  #                                          "output/calib/",
  #                                          "out_20260113_dbhXbatot_phylogeny.Rdata")),


  # COMPUTE OUTPUT VARIABLES ----

  # Create parameters table ----
  tar_target(output_params, get_output_params(models_output_list,
                                              n_analysis = 5000)),


  # ## Compute tree-level variables ---- 
  # tar_target(data_output_tree, compute_output_tree(data_calib,
  #                                                  models_setup,
  #                                                  output_params,
  #                                                  LAD_CONTROL)),
  # 
  # 
  # ## Compute stand-level variables ---- 
  # tar_target(data_output_stand, compute_output_stand(data_calib,
  #                                                    models_setup,
  #                                                    output_params,
  #                                                    LAD_CONTROL)),
  # 
  # 
  # ## Apply SamsaRalight on output stands ----
  # tar_target(data_output_light_list, compute_output_light(data_calib,
  #                                                         data_rad,
  #                                                         init_db$plots,
  #                                                         models_setup,
  #                                                         output_params,
  #                                                         LAD_CONTROL)),
  # 
  # tar_target(data_output_light, bind_output_light(data_output_light_list)),
  
  NULL
)
