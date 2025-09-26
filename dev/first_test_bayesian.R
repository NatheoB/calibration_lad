# IMPORTS ----
library(targets)
library(dplyr)
library(SamsaRaLight)
library(BayesianTools)
library(purrr)

# FUNCTIONS ----

## FUNCTION TO RUN SAMSARALIGHT ----
run_sl_standXlad <- function(data_sl, 
                             data_rad, 
                             site, 
                             rep,
                             p_lad_intercept,
                             p_lad_dbh) {
  
  # Get tree dataset and set the LAD value for all trees
  tmp_trees <- data_sl[[site]][[rep]]$trees %>% 
    dplyr::mutate(crown_lad = p_lad_intercept + dbh_cm * p_lad_dbh)
  
  # Run SamsaraLight
  tmp_out_samsalight <- 
    sl_run(tmp_trees, 
           data_rad[[site]],
           sensors = data_sl[[site]][[rep]]$sensors, 
           sensors_only = TRUE,
           latitude = 46, slope = 6, 
           aspect = 144, north_to_x_cw = 54,
           start_day = 121, end_day = 273,
           cell_size = data_sl[[site]][[rep]]$info$cell_size, 
           n_cells_x = data_sl[[site]][[rep]]$info$n_cells_x, 
           n_cells_y = data_sl[[site]][[rep]]$info$n_cells_y,
           turbid_medium = TRUE,
           trunk_interception = FALSE,
           soc = TRUE,
           height_anglemin = 15,
           direct_startoffset = 0, # =directAngleStep / 2 by default, but =0 for samsara2
           direct_anglestep = 5,
           diffuse_anglestep = 15)
  
  # Sl sensors output for the plot
  tmp_out_samsalight$sensors
  
}

## FUNCTION TO RUN SAMSARALIGHT ON ALL STANDS ----
compute_pacl_residuals <- function(data_sl,
                                   data_rad,
                                   data_sensors,
                                   exp_design,
                                   p_lad_intercept,
                                   p_lad_dbh) {
  
  # Output list
  out_residuals_list <- vector("list", nrow(exp_design))
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nrow(exp_design),
                       style = 3, width = 50, char = "=")
  
  # Run for each simulation of the experimental design
  for (i in 1:nrow(exp_design)) {
    
    # Get info about the simulation
    tmp_simu_site <- as.character(exp_design$inv_name[i])
    tmp_simu_rep <- exp_design$replicate[i]
    
    # Sl sensors output of all the plots
    tmp_out_sl <- run_sl_standXlad(data_sl, 
                                   data_rad, 
                                   tmp_simu_site, 
                                   tmp_simu_rep,
                                   p_lad_intercept,
                                   p_lad_dbh)
    
    # Compute mean residuals between all sensors
    out_residuals_list[[i]] <- dplyr::left_join(
      
      # Estimated pacl from virtual sensors
      tmp_out_sl %>%
        dplyr::select(id_sensor, pacl_sl = pacl_slope),
      
      # Measured pacl from field sensor
      data_sensors[[tmp_simu_site]] %>%
        dplyr::select(id_sensor = id, pacl_field = PACLtotal),
      
      by = "id_sensor"
    ) %>% 
      
      dplyr::mutate(residuals = pacl_sl - pacl_field) %>% 
      dplyr::select(id_sensor, pacl_sl, pacl_field, residuals)
      
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  # close(pb)
  
  dplyr::bind_rows(out_residuals_list)
}



## FUNCTION TO COMPUTE LIKELIHOOD ----
compute_likelihood <- function(p, sum = TRUE) { # sum argument necessary for computing WAIC (see BayesianTools documentation)
  
  # Parameters for the calibration run
  p_lad_intercept <- p[1] 
  p_lad_dbh <- p[2] 
  p_sigma <- p[3] # Standard deviation of residuals
  
  # Compute the residuals of the given lad parameters
  out_sl <- compute_pacl_residuals(data_sl, 
                                   data_rad, 
                                   init_db$sensors,
                                   exp_design, 
                                   p_lad_intercept,
                                   p_lad_dbh)
  residuals <- out_sl$residuals
  
  # Compute the likelihood
  # corresponds to lik 12 in Table B1 of Vrugt EMS2016 p.307
  likelihood <- dnorm(residuals, sd = p_sigma, log = TRUE) 
  
  return(ifelse(sum == TRUE,sum(likelihood),likelihood)) 
}



# SCRIPT ----

## Initial variables ----
tar_load(init_db)
tar_load(data_sl)
tar_load(data_rad)
tar_load(inv_names)
tar_load(n_replicates)

exp_design <- expand.grid(
  inv_name = inv_names,
  replicate = 1:n_replicates
)


## Get the main species to calib ----
sp_codes_occ <- purrr::map(init_db$trees, 
                           ~ table(.x$SpCode) %>% 
                             as.data.frame() %>% 
                             dplyr::rename(SpCode_SamsaraLL = Var1) %>%  
                             dplyr::mutate(SpCode_SamsaraLL = as.double(as.character(SpCode_SamsaraLL))) %>% 
                             dplyr::left_join(init_db$species[[unique(.x$plot)]], 
                                              by = "SpCode_SamsaraLL")) %>% 
  dplyr::bind_rows(.id = "site") %>% 
  dplyr::group_by(Essence_Latin) %>% 
  dplyr::summarise(Freq = sum(Freq)) %>% 
  dplyr::arrange(desc(Freq))

sp_codes_occ


## Tests ----

# ### Test the SmasaraLight function
# out_residuals <- compute_pacl_residuals(data_sl,
#                                         data_rad,
#                                         init_db$sensors,
#                                         exp_design,
#                                         p_lad_intercept = 0.7,
#                                         p_lad_dbh = 0)
# 
# ### Test the likelihood function
# out_LL <- compute_likelihood(p = c(0.7, 0, 5))


## Prior definition ----
sp_names <- 
parNames <- c("lad", "dbh", "sigma")
LB <- c(0.1, #intercept
        -0.001, #dbh effect
        0.01) #sigma
UB <- c(1.5, #lad
        0.001, #dbh effect
        20) #sigma

prior <- BayesianTools::createUniformPrior(lower = LB, upper = UB)

## Bayesian setup ----
bayesianSetup <- BayesianTools::createBayesianSetup(compute_likelihood, 
                                                    prior, 
                                                    names = parNames)

settings <- list(iterations = 15000, nCR = 3, gamma = NULL, eps = 0, e = 0.05, 
                 pCRupdate = FALSE, updateInterval = 100, burnin = 2500, thin = 1, 
                 adaptation = 0.2, parallel = NULL, Z = NULL, ZupdateFrequency = 10, 
                 pSnooker = 0.1, DEpairs = 2, consoleUpdates = 100, startValue = NULL, 
                 currentChain = 1, message = TRUE)


## Run Bayesian calibration ----
t_start <- Sys.time()
out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup, 
                              sampler = "DREAMzs", 
                              settings = settings)
t_end <- Sys.time()

computation_time <- t_end - t_start
computation_time


## Observe output ----
plot(out)
summary(out)


## Save Bayesian output ----
session <- Sys.info()
thetime <- format(Sys.time(), "%y%m%d")
save(out,prior,bayesianSetup,settings,session,computation_time,
     file = file.path(getwd(), "output", "calib", paste0("res_",thetime,".Rdata")))
