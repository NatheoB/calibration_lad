# _targets.R file
library(targets)

# library(future)
# library(future.callr)
# plan(callr)

# Source functions in R folder
lapply(grep("R$", list.files("R", recursive = TRUE), value = TRUE), function(x) source(file.path("R", x)))

# Set options (i.e. clustermq.scheduler for multiprocess computing)
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")

tar_option_set(packages = c("dplyr", "tidyr", "data.table", "vroom",
                            "SamsaRaLight"))

# install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type="source")
# install.packages("https://sourceforge.net/projects/rcapsis.capsisbridge.p/files/latest", repos = NULL,  type="source")
# devtools::install_github("ecoinfor/U.Taxonstand")
# devtools::install_github("EnquistLab/RTNRS")

# List of targets
list(
  
  # GLOBAL PARAMETERS ----
  tar_target(SEED, 5000),
  
  
  # Inputs ----
  
  ## Sensors dataset ----
  tar_target(data_sensors, expand.grid(
    x = seq(5, 95, by = 10),
    y = seq(5, 95, by = 10)
  ) %>%
    dplyr::mutate(id = row_number(),
                  h_m = 1.3)),
  
  
  ## Radiation dataset ----
  tar_target(data_rad, SamsaRaLight::data_rad_prenovel),
  
  
  ## Inventories for SamsareaLight ----
  tar_target(trees_list, list("prenovel" = SamsaRaLight::data_trees_prenovel)),
  tar_target(lads, seq(0.1, 1.5, by = 0.1)),
  tar_target(species, c("Abies alba", "Picea abies", "Fagus sylvatica")),
  
  tar_target(exp_design, 
             create_experimental_design(names(trees_list)[1], lads, species)
  ),
  
  tar_target(invs_for_samsalight, 
             create_samsalight_inventories(trees_list, exp_design)),
  
  
  # Run SamsaraLight for each LAD ----
  tar_target(out_samsalight_list, run_samsalight_expdesign(invs_for_samsalight,
                                                           data_sensors, data_rad)),
  
  
  NULL
  )
