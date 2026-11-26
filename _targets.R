# _targets.R file
library(targets)

library(future)
library(future.callr)
plan(callr)

# Install SamsaraLight
# install.packages("SamsaRaLight_1.1.0.tar.gz",
#                  repos = NULL, type = "source")

# Source functions in R folder
lapply(grep("R$", list.files("R", recursive = TRUE), value = TRUE), 
       function(x) source(file.path("R", x)))


# Set options (i.e. clustermq.scheduler for multiprocess computing)
options(tidyverse.quiet = TRUE,
        dplyr.summarise.inform = FALSE,
        clustermq.scheduler = "multiprocess")

tar_option_set(packages = c("dplyr", "tidyr", "data.table", 
                            "vroom", "purrr", "janitor",
                            "ggplot2", "ggforce", "sf", "sp",
                            "SamsaRaLight", "readxl", "nlme",
                            "BayesianTools", "extraDistr", "loo",
                            "rlang", "e1071", "mgcv", "posterior"))

# List of targets
list(
  
  # PARAMETERS ----
  tar_target(SEED, 5030),
  
  tar_target(SL_VERSION, {
    v <- packageVersion("SamsaRaLight")
    if (v < "1.1.0.9001") stop("SamsaRaLight >= 1.1.0.9001 is required (installed: ", v, ")")
    as.character(v)
  }, cue = tar_cue(mode = "always")),
  
  tar_target(LAD_CONTROL, 0.5),
  
  tar_target(N_SL_THREADS, 1), # OpenMP threads per site (benchmark 2026-11-24: workers scale better than threads)
  tar_target(N_SITE_WORKERS, 1), # Workers per chain running the sites in parallel (1 = sequential); with the SamsaRaLight geometry (2026-11-25), sites are too fast for a cluster: parallelise chains instead
  tar_target(USE_SL_GEOMETRY, TRUE), # Estimate PACL from the precomputed geometry (FALSE = run SamsaraLight at each evaluation)
  
  tar_target(N_REPS, 2), # TEST (production: 3)
  tar_target(N_ITERATIONS_PER_REP, 500), # TEST (production: 600000)
  tar_target(N_BURNING_PER_REP, 300), # TEST (production: 500000)
  tar_target(N_SAMPLES_PER_REP, 50), # TEST (production: 1000)
  tar_target(N_SAMPLES_POINTWISE_PER_REP, 50), # TEST (production: 1000)
  
  tar_target(LOCAL_BA_RADIUS, 15),
  tar_target(CELL_SIZE, 5),
  tar_target(SENSOR_TRUNK_MARGIN, 0.10), # Distance (m) from the stem surface where sensors mapped inside a trunk are moved
  
  
  # PREPARE DATABASE ----
  
  ## Format the different datasets into a calibration database ----
  tar_target(init_db_raw, format_database(
    data_irres_fp = "data/inventories/IRRES/backupIRRES_20261006_999-final.RData",
    data_cloture_fp = "data/inventories/Cloture/cloture.rdata",
    data_baileux_folderpath = "data/inventories/UCLouvain/",
    data_fontfreyde_folderpath = "data/inventories/INRAe/Fontfreyde",
    data_lorris_folderpath = "data/inventories/INRAe/Lorris",
    data_ventoux_folderpath = "data/inventories/INRAe/Ventoux",
    dico_species_fp = "data/inventories/dicoEssence.xlsx"
  )),
  
  
  ## Information about plots ----
  tar_target(plot_infos_fp, "data/inventories/plots_info.csv", format = "file"),
  tar_target(data_plots, vroom::vroom(plot_infos_fp)),

  ## Finalise the database ----
  tar_target(init_db, finalise_database(init_db_raw, data_plots)),
  
  
  # RECONSTRUCT 3D stands ----
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(data_plots)),
  tar_target(data_rad_fp, save_input_rad(data_rad, fp = "data/data_rad.csv")),
  # tar_target(data_rad, load_radiation_dataset("data/data_rad.csv")),
  
  
  ## Calibrate allometric relationships ----
  tar_target(params_alloms, fit_crown_allometries(init_db, output_folder = "output/alloms")),
  
  
  ## Create SamsaRaLight stands ----  
  tar_target(data_stands, {
    SL_VERSION # Rebuild the stands when SamsaRaLight package version changes
    create_sl_stands(init_db, 
                     data_rad,
                     params_alloms,
                     cell_size = CELL_SIZE,
                     local_ba_radius = LOCAL_BA_RADIUS,
                     sensor_trunk_margin = SENSOR_TRUNK_MARGIN,
                     seed = SEED)
  }),
  
  
  ## Precompute the SamsaRaLight geometry of each stand ----
  tar_target(sl_geometries, create_sl_geometries(data_stands, data_rad)),
  
  
  ## Plot the stands ----
  # tar_target(plots_stands_fp, plot_stands(data_stands,
  #                                         species_colors_fp = "data/inventories/species_colors.csv",
  #                                         output_plots_fp = "output/reconstructed_plots")),
  
  
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
  tar_target(id_models, c(1, 5)), # TEST (production: exp_design$id_model)


  ## Plot information in the same order as data_stands ----
  tar_target(data_plots_calib, data_plots %>%
               dplyr::slice(match(names(data_stands), site)) %>%
               dplyr::rename(name = site)),

  tar_target(site_names, names(data_stands)),


  ## Initialise the Bayesian setups ----
  tar_target(models_setup, initialise_models(exp_design,
                                             data_plots_calib,
                                             data_stands,
                                             data_rad,
                                             sensors_punobs,
                                             prior_lad = LAD_CONTROL,
                                             n_threads = N_SL_THREADS,
                                             use_sl_geometry = USE_SL_GEOMETRY,
                                             sl_geometries = sl_geometries)),


  ## Run the MCMC ----
  tar_target(reps, 1:N_REPS),

  tar_target(models_output_list, calibrate_models(models_setup,
                                                  id_model = id_models,
                                                  i_rep = reps,
                                                  n_site_workers = N_SITE_WORKERS,
                                                  logs_folder = "logs/calib"),
             pattern = cross(id_models, reps),
             iteration = "list"),

  ## Save output model ----
  tar_target(output_models_fp, save_output_models(exp_design,
                                                  models_setup,
                                                  models_output_list,
                                                  "output/calib/",
                                                  "outmods_20261126_test.Rdata")),

  
  # COMPARE MODELS ----
  
  ### Compute pointwise likelihoods ----
  tar_target(models_summary_pointwise_list, get_summary_pointwise_models(models_setup,
                                                                         models_output_list,
                                                                         N_BURNING_PER_REP,
                                                                         N_SAMPLES_POINTWISE_PER_REP,
                                                                         N_SITE_WORKERS,
                                                                         logs_folder = "logs/pointwise"),
             pattern = map(models_output_list),
             iteration = "list"),



  ### Compare models with LOO-CV and WAIC ----
  tar_target(models_comparison, compare_models(models_summary_pointwise_list)),

  ### Convergence diagnostics (Rhat, ESS) on the full post-burn-in chains ----
  tar_target(models_convergence, diagnose_convergence(models_output_list, N_BURNING_PER_REP)),


  ### Evaluate models with RMSE ----
  tar_target(models_evaluation, evaluate_models(models_summary_pointwise_list)),


  ### Save model comparison ----
  tar_target(output_comparison_fp, save_output_comparison(models_summary_pointwise_list,
                                                          models_comparison,
                                                          models_evaluation,
                                                          models_convergence,
                                                          "output/comparison/",
                                                          "outcomp_20261116_test.Rdata")),


  # # COMPUTE OUTPUT VARIABLES ----
  # tar_target(analysis_id_mods, c(1, 5)),
  # 
  # ## Create parameters table ----
  # tar_target(output_params, get_output_params(models_output_list,
  #                                             n_burning = N_BURNING_PER_REP,
  #                                             n_samples_per_chain = N_SAMPLES_PER_REP)),
  # 
  # ## Estimate LAD ----
  # tar_target(output_tree_lad, compute_output_tree_lad(data_stands,
  #                                                     models_setup,
  #                                                     output_params,
  #                                                     id_mod = analysis_id_mods,
  #                                                     site_name = site_names),
  #            pattern = cross(analysis_id_mods, site_names), iteration = "list"),
  # 
  # 
  # ## Get intercepted light with SamsaraLight ----
  # tar_target(output_tree_light, compute_output_tree_light(data_stands,
  #                                                         data_rad,
  #                                                         output_tree_lad),
  #            pattern = map(output_tree_lad), iteration = "list"),
  # 
  # 
  # ## Save output data for analysis ----
  # tar_target(output_analysis_fp, save_output_analysis(output_params,
  #                                                     output_tree_lad,
  #                                                     output_tree_light,
  #                                                     "output/analysis/",
  #                                                     "outanalysis_20261116_test.Rdata")),
  
  NULL
)
