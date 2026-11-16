fit_crown_allometries <- function(init_db, output_folder) {
  
  # Bind all trees ----
  data_trees <- init_db %>% 
    purrr::map(
      ~.x$trees %>% 
        dplyr::mutate(
          id_tree = as.character(id_tree),
          origin = .x$info$origin
        ) %>% 
        dplyr::relocate(origin)
    ) %>% 
    dplyr::bind_rows(.id = "site")
  
  
  
  # HEIGHT - DBH relationships ----
  message("Fitting tree height...")
  
  ## Prepare height dataset ----
  data_height <- data_trees %>% 
    dplyr::filter(
      !is.na(h_m),
      use_for_allom
    )
  
  ## Get species to fit ----
  sort(table(data_height$sp_allom), decreasing = TRUE)
  sp_allom_height <- unique(data_height$sp_allom)
  
  ## Fit species-specific models ----
  models_height <- fit_height_models(
    data_height, 
    sp_allom_height, 
    sp_noweight = c("Carpinus_betulus", "other_angio"),
    sp_noranef = NULL
  )
  
  ## Plot output models ----
  plot_output_height(models_height, output_folder)
  
  
  
  # CROWN BASE HEIGHT - HEIGHT relationships ----
  
  ## Prepare crown base height dataset ----
  data_crownbase <- data_trees %>% 
    dplyr::filter(
      !is.na(hbase_m),
      use_for_allom
    )
  
  ## Get species to fit ----
  sort(table(data_crownbase$sp_allom), decreasing = TRUE)
  sp_allom_crownbase <- unique(data_crownbase$sp_allom)
  
  
  ## Fit species-specific models ----
  models_crownbase <- fit_crownbase_models(
    data_crownbase, 
    sp_allom_crownbase, 
    sp_noweight = c("Abies_alba"),
    sp_noranef = c("Betula_sp")
  )
  
  
  ## Plot output models ----
  plot_output_crownbase(models_crownbase, output_folder)
  
  
  
  # CROWN RADIUS - DBH relationships ----
  
  ## Prepare crown radius dataset ----
  data_cradius <- prepare_data_cradius(data_trees)

  ## Get species to fit ----
  sort(table(data_cradius$sp_allom), decreasing = TRUE)
  sp_allom_cradius <- unique(data_cradius$sp_allom)
  
  ## Fit species-specific models ----
  models_cradius <- fit_cradius_models(
    data_cradius, 
    sp_allom_cradius
  )
  
  ## Plot output models ----
  plot_output_cradius(models_cradius, output_folder)
  
  
  
  # PARAMETERS TABLE ----
  params_alloms <- create_parameters_table(
    models_height, 
    models_crownbase, 
    models_cradius,
    output_folder
  )
  
  return(params_alloms)
}





# FUNCTIONS HEIGHT ----

fit_height_model <- function(data, sp, no_weight = FALSE, no_ranef = FALSE,
                             start_a1 = 50, start_a2 = 0.05, start_a3 = 1,
                             n_pred = 200) {
  
  # 1. Prepare data
  dat <- data %>%
    filter(
      sp_allom == sp,
      !is.na(h_m),
      !is.na(dbh_cm),
      !is.na(site),
      dbh_cm > 0
    ) %>%
    droplevels()
  
  
  # 2. Basic checks
  if (nrow(dat) == 0) {
    stop(
      paste0(
        "No valid observations found for species: ",
        sp
      )
    )
  }
  
  n_sites <- dplyr::n_distinct(dat$site)
  
  if (n_sites < 2) {
    warning(
      paste0(
        "Only ", n_sites,
        " site(s) found. The random effect of site ",
        "may not be estimable reliably."
      )
    )
  }
  
  
  # Print data
  cat("\n")
  cat("Fitting HEIGHT - DBH relationship...\n")
  cat("============================================\n")
  cat("Species:", sp, "\n")
  cat("============================================\n")
  cat("Number of trees:", nrow(dat), "\n")
  cat("Number of sites:", n_sites, "\n")
  cat("============================================\n")
  cat("DBH weight:", !no_weight, " ------ plot random effect:", !no_ranef, "\n\n")
  
  
  # 3. Fit Chapman–Richards model
  if (no_ranef) {
    
    # No random effects because it does not converge otherwise
    m0 <- gnls(
      h_m ~ 1.3 + a1 * (1 - exp(-a2 * dbh_cm))^a3,
      data = dat,
      params = list(a1 + a2 + a3 ~ 1),
      start=c(
        a1 = start_a1,
        a2 = start_a2,
        a3 = start_a3
      ),
      weights = if (no_weight) NULL else varPower(0.5,~dbh_cm)
    )
    
  } else {
    
    m0 <- nlme(
      h_m ~ 1.3 + a1 * (1 - exp(-a2 * dbh_cm))^a3,
      
      data = dat,
      
      fixed = a1 + a2 + a3 ~ 1,
      random = a1 ~ 1 | site,
      
      start = c(
        a1 = start_a1,
        a2 = start_a2,
        a3 = start_a3
      ),
      weights = if (no_weight) NULL else varPower(0.5,~dbh_cm)
    )
    
  }
  
  
  # 4. Residuals
  dat_resid <- dat %>%
    mutate(
      resid = residuals(m0)
    )
  
  mae <- mean(abs(dat_resid$resid), na.rm = TRUE)
  
  p_resid <- ggplot(
    dat_resid,
    aes(
      x = dbh_cm,
      y = resid
    )
  ) +
    geom_point(
      alpha = 0.3,
      size = 1.8
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_smooth(
      formula = 'y ~ x',
      method = "loess",
      se = TRUE
    ) +
    labs(
      x = "DBH (cm)",
      y = "Height (m) - residuals",
      title = "Residual diagnostics",
      subtitle = paste0(
        sp, " (MAE = ", round(mae, 2), " m)"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 5. Prediction grid by site
  pred_grid <- expand.grid(
    dbh_cm = seq(
      min(dat$dbh_cm, na.rm = TRUE),
      max(dat$dbh_cm, na.rm = TRUE),
      length.out = n_pred
    ),
    site = unique(dat$site)
  )
  
  pred_grid$height_pred <- predict(
    m0,
    newdata = pred_grid,
    level = 1
  )
  
  
  # 6. Population-level predictions
  pred_pop <- data.frame(
    dbh_cm = seq(
      min(dat$dbh_cm, na.rm = TRUE),
      max(dat$dbh_cm, na.rm = TRUE),
      length.out = n_pred
    )
  )
  
  pred_pop$height_pred <- predict(
    m0,
    newdata = pred_pop,
    level = 0
  )
  
  
  # 7. Prediction plot
  p_pred <- ggplot() +
    
    # Observations
    geom_point(
      data = dat,
      aes(
        x = dbh_cm,
        y = h_m,
        colour = site
      ),
      alpha = 0.35,
      size = 1.8
    ) +
    
    # Site-specific predictions
    geom_line(
      data = pred_grid,
      aes(
        x = dbh_cm,
        y = height_pred,
        colour = site,
        group = site
      ),
      linewidth = 1
    ) +
    
    # Population-level prediction
    geom_line(
      data = pred_pop,
      aes(
        x = dbh_cm,
        y = height_pred
      ),
      colour = "black",
      linewidth = 1.4,
      linetype = "dashed"
    ) +
    
    # H = 1.3 m at DBH = 0
    geom_hline(
      yintercept = 1.3,
      linetype = "dotted",
      colour = "grey40"
    ) +
    
    labs(
      x = "DBH (cm)",
      y = "Tree height (m)",
      colour = "Site",
      title = "Height – DBH relationship",
      subtitle = paste0(sp, " (", nrow(dat), " observations)")
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      legend.position = "right",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 8. Parameters
  if (no_ranef) {
    fixed_coef <- coefficients(m0)
    random_coef <- NULL
    
  } else {
    fixed_coef <- fixef(m0)
    random_coef <- ranef(m0)
  }
  
  params <- data.frame(
    site = c("POPULATION", rownames(random_coef)),
    a1 = fixed_coef[["a1"]] + c(0, random_coef[["a1"]]),
    a2 = fixed_coef[["a2"]],
    a3 = fixed_coef[["a3"]]
  )
  
  # 9. Return everything
  return(
    list(
      sp = sp,
      data = dat,
      model = m0,
      residual_data = dat_resid,
      mae = mae,
      predictions_site = pred_grid,
      predictions_population = pred_pop,
      parameters = params,
      plot_residuals = p_resid,
      plot_predictions = p_pred
    )
  )
}



fit_height_models <- function(data_height, species_names, sp_noweight, sp_noranef) {
  
  models_height <- setNames(
    vector("list", length(species_names)), 
    species_names
  )
  
  non_fitted_sp <- c()
  
  for (sp in species_names) {
    
    tmp_out <- tryCatch(
      
      fit_height_model(
        data_height, sp, 
        no_weight = sp %in% sp_noweight,
        no_ranef = sp %in% sp_noranef
      ),
      
      error = function(e) {
        message(sp, " ---- ", e)
        NULL
      }
    )
    
    if (is.null(tmp_out)) {
      non_fitted_sp <- c(non_fitted_sp, sp)
    } else {
      models_height[[sp]] <- tmp_out
    }
    
  }
  
  message("Species that have not converged: ", paste(non_fitted_sp, collapse = ", "))
  
  models_height
}



plot_output_height <- function(models_height, output_folder) {
  
  message("Plotting height allometries....")
  
  
  # 1. Create height output folder
  height_folder <- file.path(
    output_folder,
    "height"
  )
  
  dir.create(
    height_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # 2. Species-specific plots
  for (sp in names(models_height)) {
    
    if (is.null(models_height[[sp]])) next
    
    # Residual plot
    if (!is.null(models_height[[sp]]$plot_residuals)) {
      
      ggplot2::ggsave(
        filename = file.path(
          height_folder,
          paste0(
            "height_residuals_",
            sp,
            ".png"
          )
        ),
        plot = models_height[[sp]]$plot_residuals,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
    
    
    # Prediction plot
    if (!is.null(models_height[[sp]]$plot_predictions)) {
      
      ggplot2::ggsave(
        filename = file.path(
          height_folder,
          paste0(
            "height_predictions_",
            sp,
            ".png"
          )
        ),
        plot = models_height[[sp]]$plot_predictions,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
  
  
  # 3. Population comparison plot
  predictions_population <- models_height %>%
    purrr::map(
      ~ .x$predictions_population
    ) %>%
    dplyr::bind_rows(
      .id = "species"
    )
  
  
  plot_comparison <- predictions_population %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = dbh_cm,
        y = height_pred,
        color = species
      )
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.1
    ) +
    
    ggplot2::geom_hline(
      yintercept = 1.3,
      linetype = "dotted",
      colour = "grey40"
    ) +
    
    ggplot2::labs(
      x = "DBH (cm)",
      y = "Tree height (m)",
      colour = "Species",
      title = "Height–DBH relationship",
      subtitle = "Within each species' observed DBH range"
    ) +
    
    ggplot2::theme_classic(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 4. Save population comparison
  ggplot2::ggsave(
    filename = file.path(
      height_folder,
      "height_predicitons.png"
    ),
    plot = plot_comparison,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  
  invisible(plot_comparison)
}



# FUNCTIONS CROWN BASE HEIGHT ----

fit_crownbase_model <- function(data, sp, 
                                no_weight = FALSE, no_ranef = FALSE,
                                start_a = NULL, start_b = NULL,
                                n_pred = 200) {
  
  # 1. Prepare data
  dat <- data %>%
    filter(
      sp_allom == sp,
      !is.na(hbase_m),
      !is.na(h_m),
      !is.na(site),
      !is.na(id_tree),
      hbase_m > 0,
      h_m > 0
    ) %>%
    droplevels() %>% 
    dplyr::mutate(
      id_tree = as.factor(id_tree)
    )
  
  
  # 2. Basic checks
  if (nrow(dat) == 0) {
    stop(
      paste0(
        "No valid observations found for species: ",
        sp
      )
    )
  }
  
  n_sites <- n_distinct(dat$site)
  n_trees <- n_distinct(
    interaction(dat$site, dat$id_tree, drop = TRUE)
  )
  
  
  if (n_sites < 2) {
    warning(
      paste0(
        "Only ", n_sites,
        " site(s) found. The plot-level random effect ",
        "may not be estimable reliably."
      )
    )
  }
  
  
  # Print data
  cat("\n")
  cat("Fitting CROWN BASE HEIGHT - HEIGHT relationship...\n")
  cat("============================================\n")
  cat("Species:", sp, "\n")
  cat("============================================\n")
  cat("Number of trees:", n_trees, "\n")
  cat("Number of sites:", n_sites, "\n")
  cat("============================================\n")
  cat("Height weight:", !no_weight, " ------ plot random effect:", !no_ranef, "\n\n")
  
  
  # 3. Fit power-law model
  #
  # Crown base height = a * height^b
  lm_start <- lm(
    log(hbase_m) ~ log(h_m),
    data = dat
  )
  
  if (is.null(start_a)) start_a <- exp(coef(lm_start)[1])
  if (is.null(start_b)) start_b <- coef(lm_start)[2]
  
  if (no_ranef) {
    
    # No random effects because it does not converge otherwise
    m0 <- gnls(
      hbase_m ~ a * h_m ^ b,
      data = dat,
      params = list(a + b ~ 1),
      start = c(
        a = start_a, 
        b = start_b
      ),
      weights = if (no_weight) NULL else varPower(0.5,~h_m)
    )
    
  } else {
    
    m0 <- nlme(
      hbase_m ~ a * h_m ^ b,
      data = dat,
      fixed = a + b ~ 1,
      random = a ~ 1 | site,
      start=c(
        a = start_a,
        b = start_b
      ),
      weights = if (no_weight) NULL else varPower(0.5,~h_m)
    )
    
  }

  
  
  # 4. Residuals
  dat_resid <- dat %>%
    mutate(
      log_crown = log(h_m),
      resid = residuals(m0)
    )
  
  mae <- mean(abs(dat_resid$resid), na.rm = TRUE)
  
  p_resid <- ggplot(
    dat_resid,
    aes(
      x = h_m,
      y = resid
    )
  ) +
    geom_point(
      alpha = 0.3,
      size = 1.8
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_smooth(
      formula = 'y ~ x',
      method = "loess",
      se = TRUE
    ) +
    labs(
      x = "Tree height (m)",
      y = "Crown base height (m) - residuals",
      title = "Crown base height model — residual diagnostics",
      subtitle = paste0(
        sp, " (MAE = ", round(mae, 2), " m)"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 5. Parameters
  if (no_ranef) {
    fixed_coef <- coefficients(m0)
    random_coef <- NULL
    
  } else {
    fixed_coef <- fixef(m0)
    random_coef <- ranef(m0)
  }
  
  params <- data.frame(
    site = c("POPULATION", rownames(random_coef)),
    a = fixed_coef[["a"]] + c(0, random_coef[["a"]]),
    b = fixed_coef[["b"]]
  )
  
  
  # 6. Prediction grid
  pred_grid <- expand.grid(
    
    h_m = seq(min(dat$h_m, na.rm = TRUE),
                   max(dat$h_m, na.rm = TRUE),
                   length.out = 200),
    
    site = params$site
  ) %>% 
    dplyr::left_join(
      params, by = "site"
    ) %>% 
    dplyr::mutate(
      crownbase = a * h_m^b
    )
  
  pred_plot <- pred_grid %>% dplyr::filter(site != "POPULATION")
  pred_pop <- pred_grid %>% dplyr::filter(site == "POPULATION")
  
  
  # 7. Prediction plot
  p_pred <- ggplot() +
    
    # Observed measurements
    geom_point(
      data = dat,
      aes(
        x = h_m,
        y = hbase_m,
        colour = site,
        group = site
      ),
      alpha = 0.25,
      size = 1.5
    ) +
    
    # Plot-specific curves
    geom_line(
      data = pred_plot,
      aes(
        x = h_m,
        y = crownbase,
        colour = site,
        group = site
      ),
      linewidth = 1
    ) +
    
    # Population curve
    geom_line(
      data = pred_pop,
      aes(
        x = h_m,
        y = crownbase
      ),
      colour = "black",
      linewidth = 1.3,
      linetype = "dashed"
    ) +
    
    labs(
      x = "Tree height (m)",
      y = "Crown base height (m)",
      colour = "Plot",
      title = "Crown base height - Tree height relationships",,
      subtitle = paste0(sp, " (", nrow(dat), " observations)")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  # 9. Return results
  return(
    list(
      sp = sp,
      data = dat,
      model = m0,
      parameters = params,
      residual_data = dat_resid,
      mae = mae,
      predictions_plot = pred_plot,
      predictions_population = pred_pop,
      plot_predictions = p_pred,
      plot_residuals = p_resid
    )
  )
}


fit_crownbase_models <- function(data_crownbase, species_names, 
                                 sp_noweight = NULL, sp_noranef = NULL) {
  
  models_crownbase <- setNames(
    vector("list", length(species_names)), 
    species_names
  )
  
  non_fitted_sp <- c()
  
  for (sp in species_names) {
    
    tmp_out <- tryCatch(
      
      fit_crownbase_model(
        data_crownbase, sp, 
        no_weight = sp %in% sp_noweight,
        no_ranef = sp %in% sp_noranef
      ),
      
      error = function(e) {
        message(sp, " ---- ", e)
        NULL
      }
    )
    
    if (is.null(tmp_out)) {
      non_fitted_sp <- c(non_fitted_sp, sp)
    } else {
      models_crownbase[[sp]] <- tmp_out
    }
    
  }
  
  message("Species that have not converged: ", paste(non_fitted_sp, collapse = ", "))
  
  models_crownbase
}


plot_output_crownbase <- function(models_crownbase, output_folder) {
  
  message("Plotting crown base height allometries....")
  
  # 1. Create crown base output folder
  crownbase_folder <- file.path(
    output_folder,
    "crownbase"
  )
  
  dir.create(
    crownbase_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # 2. Species-specific plots
  for (sp in names(models_crownbase)) {
    
    if (is.null(models_crownbase[[sp]])) next
    
    # Residual plot
    if (!is.null(models_crownbase[[sp]]$plot_residuals)) {
      
      ggplot2::ggsave(
        filename = file.path(
          crownbase_folder,
          paste0(
            "crownbase_residuals_",
            sp,
            ".png"
          )
        ),
        plot = models_crownbase[[sp]]$plot_residuals,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
    
    
    # Prediction plot
    if (!is.null(models_crownbase[[sp]]$plot_predictions)) {
      
      ggplot2::ggsave(
        filename = file.path(
          crownbase_folder,
          paste0(
            "crownbase_predictions_",
            sp,
            ".png"
          )
        ),
        plot = models_crownbase[[sp]]$plot_predictions,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
  
  
  # 3. Population comparison plot
  predictions_population <- models_crownbase %>%
    purrr::map(
      ~ .x$predictions_population
    ) %>%
    dplyr::bind_rows(
      .id = "species"
    )
  
  
  plot_comparison <- predictions_population %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = h_m,
        y = crownbase,
        color = species
      )
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.1
    ) +
    
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      colour = "grey40"
    ) +
    
    ggplot2::labs(
      x = "Tree height (m)",
      y = "Crown base height (m)",
      colour = "Species",
      title = "Crown base – Total height relationship",
      subtitle = "Within each species' observed range of height"
    ) +
    
    ggplot2::theme_classic(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 4. Save population comparison
  ggplot2::ggsave(
    filename = file.path(
      crownbase_folder,
      "crownbase_height_population.png"
    ),
    plot = plot_comparison,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  
  invisible(plot_comparison)
}



# FUNCTIONS CROWN RADIUS ----

prepare_data_cradius <- function(data_trees) {
  
  data_cradius <- data_trees %>%
    tidyr::pivot_longer(
      cols = matches("^(rn|re|rs|rw)_m$"),
      names_to = "radius_type",
      values_to = "radius_m"
    ) %>% 
    dplyr::filter(
      !is.na(radius_m),
      use_for_allom
    )
  
  return(data_cradius)
}


fit_cradius_model <- function(data, sp,
                              start_a = NULL, start_b = NULL,
                              n_pred = 200) {
  
  # 1. Prepare data
  dat <- data %>%
    filter(
      sp_allom == sp,
      !is.na(radius_m),
      !is.na(dbh_cm),
      !is.na(site),
      !is.na(id_tree),
      radius_m > 0,
      dbh_cm > 0
    ) %>%
    droplevels() %>% 
    dplyr::mutate(
      id_tree = as.factor(id_tree)
    )
  
  
  # 2. Basic checks
  if (nrow(dat) == 0) {
    stop(
      paste0(
        "No valid observations found for species: ",
        sp
      )
    )
  }
  
  n_sites <- n_distinct(dat$site)
  n_trees <- n_distinct(
    interaction(dat$site, dat$id_tree, drop = TRUE)
  )
  
  
  if (n_sites < 2) {
    warning(
      paste0(
        "Only ", n_sites,
        " site(s) found. The plot-level random effect ",
        "may not be estimable reliably."
      )
    )
  }
  
  # Print data
  cat("\n")
  cat("Fitting CROWN RADIUS - DBH relationship...\n")
  cat("============================================\n")
  cat("Species:", sp, "\n")
  cat("============================================\n")
  cat("Number of observations:", nrow(dat), "\n")
  cat("Number of trees:", n_trees, "\n")
  cat("Number of sites:", n_sites, "\n\n")
  
  
  
  # 3. Fit power-law model
  #
  # Crown radius = a * DBH^b
  lm_start <- lm(
    log(radius_m) ~ log(dbh_cm),
    data = dat
  )
  
  if (is.null(start_a)) start_a <- exp(coef(lm_start)[1])
  if (is.null(start_b)) start_b <- coef(lm_start)[2]
  
  m0 <- nlme(
    radius_m ~ a * dbh_cm ^ b,
    data = dat,
    fixed = a + b ~ 1,
    random = a ~ 1 | site,
    start=c(
      a = start_a,
      b = start_b
    ),
    weights = varPower(0.5,~dbh_cm)
  )
  
  
  
  # 4. Residuals
  dat_resid <- dat %>%
    mutate(
      log_crown = log(radius_m),
      resid = residuals(m0)
    )
  
  
  mae <- mean(abs(dat_resid$resid), na.rm = TRUE)
  
  p_resid <- ggplot(
    dat_resid,
    aes(
      x = dbh_cm,
      y = resid
    )
  ) +
    geom_point(
      alpha = 0.3,
      size = 1.8
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_smooth(
      formula = 'y ~ x',
      method = "loess",
      se = TRUE
    ) +
    labs(
      x = "DBH (cm)",
      y = "Crown radius (m) - residuals",
      title = "Crown radius model — residual diagnostics",
      subtitle = paste0(
        sp, " (MAE = ", round(mae, 2), " m)"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 5. Parameters
  fixed_coef <- fixef(m0)
  random_coef <- ranef(m0)
  
  params <- data.frame(
    site = c("POPULATION", rownames(random_coef)),
    a = fixed_coef[["a"]] + c(0, random_coef[["a"]]),
    b = fixed_coef[["b"]]
  )
  
  
  # 6. Prediction grid
  pred_grid <- expand.grid(
    
    dbh_cm = seq(min(dat$dbh_cm, na.rm = TRUE),
                 max(dat$dbh_cm, na.rm = TRUE),
                 length.out = 200),
    
    site = params$site
  ) %>% 
    dplyr::left_join(
      params, by = "site"
    ) %>% 
    dplyr::mutate(
      radius = a * dbh_cm^b
    )
  
  pred_plot <- pred_grid %>% dplyr::filter(site != "POPULATION")
  pred_pop <- pred_grid %>% dplyr::filter(site == "POPULATION")
  
  
  
  # 7. Prediction plot
  p_pred <- ggplot() +
    
    # Observed measurements
    geom_point(
      data = dat,
      aes(
        x = dbh_cm,
        y = radius_m,
        colour = site
      ),
      alpha = 0.25,
      size = 1.5
    ) +
    
    # Plot-specific curves
    geom_line(
      data = pred_plot,
      aes(
        x = dbh_cm,
        y = radius,
        colour = site,
        group = site
      ),
      linewidth = 1
    ) +
    
    # Population curve
    geom_line(
      data = pred_pop,
      aes(
        x = dbh_cm,
        y = radius
      ),
      colour = "black",
      linewidth = 1.3,
      linetype = "dashed"
    ) +
    
    labs(
      x = "DBH (cm)",
      y = "Crown radius (m)",
      colour = "Plot",
      title = "Crown radius - DBH relationships",,
      subtitle = paste0(sp, " (", nrow(dat), " observations)")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  # 9. Return results
  return(
    list(
      sp = sp,
      data = dat,
      model = m0,
      parameters = params,
      residual_data = dat_resid,
      mae = mae,
      predictions_plot = pred_plot,
      predictions_population = pred_pop,
      plot_predictions = p_pred,
      plot_residuals = p_resid
    )
  )
}



fit_cradius_models <- function(data_cradius, species_names) {
  
  models_cradius <- setNames(
    vector("list", length(species_names)), 
    species_names
  )
  
  non_fitted_sp <- c()
  
  for (sp in species_names) {
    
    tmp_out <- tryCatch(
      withCallingHandlers(
        fit_cradius_model(data_cradius, sp),
        warning = function(w) {
          message(sp, " ---- WARNING: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message(sp, " ---- ERROR: ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(tmp_out)) {
      non_fitted_sp <- c(non_fitted_sp, sp)
    } else {
      models_cradius[[sp]] <- tmp_out
    }
    
  }
  
  message("Species that have not converged: ", paste(non_fitted_sp, collapse = ", "))
  
  models_cradius
}



plot_output_cradius <- function(models_cradius, output_folder) {
  
  message("Plotting crown radius allometries....")
  
  # 1. Create crown radius output folder
  crownradius_folder <- file.path(
    output_folder,
    "crownradius"
  )
  
  dir.create(
    crownradius_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # 2. Species-specific plots
  for (sp in names(models_cradius)) {
    
    if (is.null(models_cradius[[sp]])) next
    
    # Residual plot
    if (!is.null(models_cradius[[sp]]$plot_residuals)) {
      
      ggplot2::ggsave(
        filename = file.path(
          crownradius_folder,
          paste0(
            "crownradius_residuals_",
            sp,
            ".png"
          )
        ),
        plot = models_cradius[[sp]]$plot_residuals,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
    
    
    # Prediction plot
    if (!is.null(models_cradius[[sp]]$plot_predictions)) {
      
      ggplot2::ggsave(
        filename = file.path(
          crownradius_folder,
          paste0(
            "crownradius_predictions_",
            sp,
            ".png"
          )
        ),
        plot = models_cradius[[sp]]$plot_predictions,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
  
  
  # 3. Population comparison plot
  predictions_population <- models_cradius %>%
    purrr::map(
      ~ .x$predictions_population
    ) %>%
    dplyr::bind_rows(
      .id = "species"
    )
  
  
  plot_comparison <- predictions_population %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = dbh_cm,
        y = radius,
        color = species
      )
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.1
    ) +
    
    ggplot2::labs(
      x = "DBH (cm)",
      y = "Crown radius (m)",
      colour = "Species",
      title = "Crown radius – DBH relationship",
      subtitle = "Within each species' observed range of DBH"
    ) +
    
    ggplot2::theme_classic(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey40",
        hjust = 0.5
      )
    )
  
  
  # 4. Save population comparison
  ggplot2::ggsave(
    filename = file.path(
      crownradius_folder,
      "crownradius_dbh_population.png"
    ),
    plot = plot_comparison,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  
  invisible(plot_comparison)
}



# FUNCTION PARAMETERS ----

create_parameters_table <- function(models_height, models_crownbase, models_cradius,
                                    output_folder) {
  
  params_alloms <- list(
    
    "height" = models_height %>% 
      purrr::map(~.x$parameters) %>% 
      dplyr::bind_rows(.id = "species"),
    
    "crownbase" = models_crownbase %>% 
      purrr::map(~.x$parameters) %>% 
      dplyr::bind_rows(.id = "species"),
    
    "cradius" = models_cradius %>% 
      purrr::map(~.x$parameters) %>% 
      dplyr::bind_rows(.id = "species")
  )
  
  write.table(
    params_alloms$height, file.path(output_folder, "params_height.csv"),
    sep = ";", dec = ".", row.names = F
  )
  
  write.table(
    params_alloms$crownbase, file.path(output_folder, "params_crownbase.csv"),
    sep = ";", dec = ".", row.names = F
  )
  
  write.table(
    params_alloms$cradius, file.path(output_folder, "params_cradius.csv"),
    sep = ";", dec = ".", row.names = F
  )
  
  return(params_alloms)
}

