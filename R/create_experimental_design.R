create_experimental_design <- function() {
  
  # 1. Control model, with and without sensors filtering ----
  # log(LAD) ~ 1 + (1 | site/origin)
  exp1 <- expand.grid(
    
    filter_sensors = c(FALSE, TRUE),
    
    site_rd_effect = TRUE,
    origin_rd_effect = TRUE,
    
    intercept_per_sp = FALSE,
    intercept_sp_pooling = FALSE,
    
    dbh_effect = FALSE,
    dbh_per_sp = FALSE,
    dbh_sp_pooling = FALSE,
    
    consider_covariance = FALSE,
    
    compet_effect = FALSE

  )  %>% 
    dplyr::mutate(
      id_exp = 1,
      id_mod = row_number()
    )
  
  
  
  # 2. Non species-specific models, adding single dbh effect ----
  # log(LAD) ~ 1 + dbh + (1 | site/origin)
  exp2 <- expand.grid(
    
    filter_sensors = TRUE,
    
    site_rd_effect = TRUE,
    origin_rd_effect = TRUE,
    
    intercept_per_sp = FALSE,
    intercept_sp_pooling = FALSE,
    
    dbh_effect = TRUE,
    dbh_per_sp = FALSE,
    dbh_sp_pooling = FALSE,
    
    consider_covariance = FALSE,
    
    compet_effect = FALSE
    
  )  %>% 
    dplyr::mutate(
      id_exp = 2,
      id_mod = row_number()
    )
  
  
  # 3. Species-specific intercept with partial pooling, with/without mean dbh ----
  # log(LAD) ~ sp + (1 | site/origin)
  # log(LAD) ~ sp + dbh + (1 | site/origin)
  exp3 <- expand.grid(
    
    filter_sensors = TRUE,
    
    site_rd_effect = TRUE,
    origin_rd_effect = TRUE,
    
    intercept_per_sp = TRUE,
    intercept_sp_pooling = TRUE,
    
    dbh_effect = c(FALSE, TRUE),
    dbh_per_sp = FALSE,
    dbh_sp_pooling = FALSE,
    
    consider_covariance = FALSE,
    
    compet_effect = FALSE
    
  )  %>% 
    dplyr::mutate(
      id_exp = 3,
      id_mod = row_number()
    )
  
  
  # 4. Species-specific intercept and dbh effect, considering covariance or not ----
  # log(LAD) ~ sp + dbh + sp:dbh + (1 | site/origin)
  exp4 <- expand.grid(
    
    filter_sensors = TRUE,
    
    site_rd_effect = TRUE,
    origin_rd_effect = TRUE,
    
    intercept_per_sp = TRUE,
    intercept_sp_pooling = TRUE,
    
    dbh_effect = TRUE,
    dbh_per_sp = TRUE,
    dbh_sp_pooling = TRUE,
    
    consider_covariance = c(FALSE, TRUE),
    
    compet_effect = FALSE
    
  )  %>% 
    dplyr::mutate(
      id_exp = 4,
      id_mod = row_number()
    )
  
  
  # FINAL EXPERIMENTAL DESIGN ----
  exp1 %>% 
    dplyr::bind_rows(exp2) %>% 
    dplyr::bind_rows(exp3) %>% 
    dplyr::bind_rows(exp4) %>% 
    
      dplyr::mutate(
        id_simu = paste0(id_exp, "_", id_mod),
        n_chains = 3,
        n_iterations = 100,
        n_analysis = 10
      ) %>%
      dplyr::relocate(id_simu, id_exp, n_chains, n_iterations, n_analysis, id_mod)
  
}
