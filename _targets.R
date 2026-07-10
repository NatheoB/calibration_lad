# _targets.R file
library(targets)

library(future)
library(future.callr)
plan(callr)

# Install SamsaraLight
# install.packages("SamsaRaLight_1.0.tar.gz",
#                  repos = NULL, type = "source")

# Source functions in R folder
lapply(grep("R$", list.files("R", recursive = TRUE), value = TRUE), 
       function(x) source(file.path("R", x)))


# Set options (i.e. clustermq.scheduler for multiprocess computing)
options(tidyverse.quiet = TRUE,
        dplyr.summarise.inform = FALSE,
        clustermq.scheduler = "multiprocess")

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
  
  tar_target(N_SL_THREADS, 8),
  
  tar_target(N_REPS, 3),
  tar_target(N_ITERATIONS_PER_REP, 600000),
  tar_target(N_BURNING_PER_REP, 500000),
  tar_target(N_SAMPLES_PER_REP, 1000),
  
  tar_target(LOCAL_BA_RADIUS, 15),
  tar_target(CELL_SIZE, 5),
  
  
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
  tar_target(data_rad_fp, save_input_rad(data_rad, fp = "data/data_rad.csv")),
  
  # tar_target(data_rad, load_radiation_dataset("data/data_rad.csv")),
  
  
  ## Create SamsaRaLight stands from tree inventories ----
  tar_target(data_stands, create_sl_stands(init_db,
                                           data_rad,
                                           cell_size = CELL_SIZE,
                                           local_ba_radius = LOCAL_BA_RADIUS,
                                           seed = SEED)),

  tar_target(site_names, names(data_stands)),
  
  
  ## Plot the calibration data ----
  # tar_target(plots_stands_fp, plot_stands(data_stands,
  #                                         data_species = init_db$species,
  #                                         output_plots_fp = "output/initial_sites")),
  
  
  
  # PRELIMINARY ANALYSIS ----
  # simple minimization of residuals
  # tar_target(lads_method1, seq(0.01, 5, by = 0.01)),
  # tar_target(lad_convergence_threshold, 3),
  # 
  # 
  # ## Estimate PACL given site LAD ----
  # tar_target(output_pacl_method1, get_sensors_pacl_sitespecificLAD(lads_method1,
  #                                                                  data_stands,
  #                                                                  data_rad)),
  # 
  # ## Optimize LAD per site ----
  # tar_target(output_lad_method1, fit_lad_method1(output_pacl_method1,
  #                                                lad_convergence_threshold,
  #                                                data_stands,
  #                                                "output/residuals_sensors")),


  # BAYESIAN CALIBRATION ----

  ## Filter sensors ----
  tar_target(sensors_punobs, filter_sensors(data_stands, 
                                            data_rad,
                                            min_pacl = 0.05,
                                            max_punobs = 0.75)),
  
  ## Create the experimental design ----
  tar_target(exp_design, create_experimental_design(N_ITERATIONS_PER_REP)),
  tar_target(id_models, exp_design$id_model),


  ## Initialise the Bayesian setups ----
  tar_target(models_setup, initialise_models(exp_design,
                                             init_db$plots,
                                             data_stands,
                                             data_rad,
                                             sensors_punobs,
                                             prior_lad = LAD_CONTROL,
                                             n_threads = N_SL_THREADS)),


  ## Run the MCMC ----
  tar_target(reps, 1:N_REPS),
  
  tar_target(models_output_list, calibrate_models(models_setup,
                                                  id_model = id_models,
                                                  i_rep = reps,
                                                  logs_folder = "logs/calib"),
             pattern = cross(id_models, reps),
             iteration = "list"),
  
  ## Save output model ----
  tar_target(output_models_fp, save_output_models(exp_design,
                                                  models_setup,
                                                  models_output_list,
                                                  "output/calib/",
                                                  "outmods_20260710_final.Rdata")),
  


  # COMPARE MODELS ----

  ### Compute pointwise likelihoods ----
  tar_target(models_summary_pointwise_list, get_summary_pointwise_models(models_setup,
                                                                         models_output_list,
                                                                         N_BURNING_PER_REP,
                                                                         logs_folder = "logs/pointwise"),
             pattern = map(models_output_list),
             iteration = "list"),



  ### Compare models with LOO-CV and WAIC ----
  tar_target(models_comparison, compare_models(models_summary_pointwise_list)),


  ### Evaluate models with RMSE ----
  tar_target(models_evaluation, evaluate_models(models_summary_pointwise_list)),

  
  ### Save model comparison ----
  tar_target(output_comparison_fp, save_output_comparison(models_summary_pointwise_list,
                                                          models_comparison,
                                                          models_evaluation,
                                                          "output/comparison/",
                                                          "outcomp_20260710_final.Rdata")),
  

  # COMPUTE OUTPUT VARIABLES ----
  tar_target(analysis_id_mods, c(2, 6)),
  
  ## Create parameters table ----
  tar_target(output_params, get_output_params(models_output_list,
                                              n_burning = N_BURNING_PER_REP,
                                              n_samples_per_chain = N_SAMPLES_PER_REP)),
  
  ## Estimate LAD ----
  tar_target(output_tree_lad, compute_output_tree_lad(data_stands,
                                                      models_setup,
                                                      output_params,
                                                      id_mod = analysis_id_mods,
                                                      site_name = site_names),
             pattern = cross(analysis_id_mods, site_names), iteration = "list"),
  

  ## Get intercepted light with SamsaraLight ----
  tar_target(output_tree_light, compute_output_tree_light(data_stands,
                                                          data_rad,
                                                          output_tree_lad),
             pattern = map(output_tree_lad), iteration = "list"),

  
  ## Save output data for analysis ----
  tar_target(output_analysis_fp, save_output_analysis(output_params,
                                                      output_tree_lad,
                                                      output_tree_light,
                                                      "output/analysis/",
                                                      "outanalysis_20260710_final.Rdata")),
  
  NULL
)
