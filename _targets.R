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
                            "BayesianTools", "extraDistr"))

# List of targets
list(
  
  # PARAMETERS ----
  
  ## Global parameters ----
  tar_target(SEED, 5030),
  
  ## SL parameters ----

  
  ## Experimental parameters ----
  tar_target(n_replicates, 1),
  
  ## Calibration parameters ----
  tar_target(n_chains, 3),
  tar_target(n_iterations, 1000),
  tar_target(n_burning, 0),

  
  
  
  # PREPARE VIRTUAL PLOTS ----
  
  ## Load and clean the initial database ----
  tar_target(init_db_fp, "data/dataBase.RData", format = "file"),
  tar_target(init_db, load_database_from_rdata(init_db_fp)),
  
  
  ## Get some informations about the database ----
  tar_target(site_names, unique(init_db$plots$name)),
  
  tar_target(sp_calib_occ, get_occurences_species2calib(init_db)),
  tar_target(species2calib, as.character(unique(sp_calib_occ$species))),
  
  
  ## Create virtual plots from tree inventories ----
  tar_target(data_sl, create_samsaralight_stands(init_db, 
                                                 n_replicates,
                                                 "output/sites",
                                                 SEED)),
  
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(init_db$plots)),
  
  
  
  
  
  # CALIBRATE THE LAD ----
  
  ## Create the experimental design ----
  tar_target(exp_design, create_experimental_design(site_names, 
                                                    n_replicates)),
  
  ## Run the Bayesian calibration with model inversion ----
  tar_target(out_calib, calibrate_lad(init_db,
                                      data_sl,
                                      data_rad,
                                      exp_design,
                                      species2calib,
                                      site_names,
                                      n_chains,
                                      n_iterations,
                                      n_burning,
                                      output_folder = file.path(getwd(), "output/calib"))),
  
  
  NULL
  )
