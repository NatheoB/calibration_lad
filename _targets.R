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
                            "SamsaRaLight"))

# List of targets
list(
  
  # PARAMETERS ----
  
  ## Global parameters ----
  tar_target(SEED, 5030),
  
  ## SL parameters ----
  tar_target(cell_size, 10),
  
  ## Experimental parameters ----
  tar_target(n_replicates, 1),
  tar_target(lad_values, seq(0.01, 2, by = 0.01)),
  
  
  
  
  # PREPARE VIRTUAL PLOTS ----
  
  ## Load and clean the initial database ----
  tar_target(init_db_fp, "data/dataBase.RData", format = "file"),
  tar_target(init_db, load_database_from_rdata(init_db_fp)),
  
  
  ## Get some informations about the database ----
  tar_target(inv_names, unique(init_db$plots$name)),
  
  
  ## Create virtual plots from tree inventories ----
  tar_target(data_sl, create_samsaralight_stands(init_db, 
                                                 cell_size, 
                                                 n_replicates,
                                                 "output/sites",
                                                 SEED)),
  
  
  ## Catch monthly radiation data from PVGIS ----
  tar_target(data_rad, get_radiation_dataset(init_db$plots)),
  
  
  
  
  
  # ESTIMATE LIGHT ON VIRTUAL SENSORS ----
  
  ## Create the experimental design ----
  tar_target(exp_design, create_experimental_design(inv_names, 
                                                    n_replicates, 
                                                    lad_values)),
  
  ## Run SamsaraLight ----
  ## On all sites, replicated as a plot, and setting a given mean LAD value
  # tar_target(out_sl, run_samsalight_expdesign(data_sl,
  #                                             data_rad,
  #                                             exp_design)),
  # 
  # 
  # 
  # 
  # 
  # # CALIBRATE LAD ----
  # 
  # ## Compute residuals ----
  # tar_target(out_residuals, compute_residuals(init_db$sensors,
  #                                             out_sl)),
  
  
  NULL
  )
