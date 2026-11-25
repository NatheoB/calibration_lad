create_sl_stands <- function(init_db,
                             data_rad,
                             params_alloms,
                             cell_size,
                             local_ba_radius,
                             seed) {
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Output list
  out_stands <- setNames(
    vector("list", length = length(init_db)),
    names(init_db)
  )
  

  # For each site 
  for (site in names(out_stands)) {
    
    message("---- Creating SamsaRaLight stand for ", site, "...")
    
    tmp_data <- init_db[[site]]
    
    # If the stand belong to IRRES2018, correct for crown measurements
    # Because PACL measurement have been taken in 2018 and crown in 2015/2016
    # Thus apply derivative of allometries to DBH increment for estimating crown measures 2 years after
    if (tmp_data$info$origin == "irres2018") {
      
      tmp_data$trees <- correct_IRRES2018_crowns(tmp_data$trees, params_alloms, site)
    }

    # Apply allometries
    tmp_data$trees <- estimate_missing_crowns(tmp_data$trees, params_alloms, site)
    
    # Set crown type and hmax
    tmp_data$trees <- set_crown_shape(tmp_data$trees, params_alloms)
    
    # Compute tree volume
    tmp_data$trees <- compute_volume_trees(tmp_data$trees)

    
    # Format and check trees inventory
    tmp_data$trees <- tmp_data$trees %>%
      dplyr::mutate(
        crown_lad = 0.5 # A default value for SL inventory
      ) %>%
      dplyr::select(
        id_tree, year,
        species = sp_calib, sp_in_inv, sp_allom,
        shadetol, leaf_group, crown_type,
        x, y, dbh_cm,
        h_m, hbase_m, hmax_m,
        rn_m, re_m, rs_m, rw_m,
        volume_m3,
        height_measured, crownbase_measured,
        hmax_measured, radii_all_measured,
        volume_measured,
        crown_lad
      )

    SamsaRaLight::check_inventory(tmp_data$trees)


    # Check plot extent
    if (!is.null(tmp_data$polygon)) {
      SamsaRaLight::check_polygon(tmp_data$polygon, tmp_data$trees)
    } else {
      message("No given core polygon: it will be automatically estimated...")
    }


    # Check sensors
    SamsaRaLight::check_sensors(tmp_data$sensors)


    # Create a square plot with the inventory in the center
    tmp_stand <- SamsaRaLight::create_sl_stand(
      trees_inv = tmp_data$trees,
      cell_size = cell_size,
      latitude = tmp_data$info$latitude,
      slope = tmp_data$info$slope,
      aspect = tmp_data$info$aspect,
      north2x = tmp_data$info$northToX,
      sensors = tmp_data$sensors,
      core_polygon_df = tmp_data$polygon,
      modify_polygon = ifelse(tmp_data$info$rect_zone, "rect", "none"),
      fill_around = TRUE
    )


    # Compute basal area competition variables
    tmp_stand$trees <- compute_competition_trees(tmp_stand$trees,
                                                 local_ba_radius,
                                                 tmp_stand$transform$new_area_ha,
                                                 tmp_stand$geometry$n_cells_x * cell_size,
                                                 tmp_stand$geometry$n_cells_y * cell_size)

    
    # Compute ray competition index RCI for each tree and percentage of unobstructed rays for each sensor
    tmp_stand <- compute_light_variables(tmp_stand, data_rad[[site]])
    
    # Add origin information in the output stand object
    tmp_stand$trees$origin <- tmp_data$info$origin
    
    # Add the virtual plot to the list
    out_stands[[site]] <- tmp_stand
  }
  
  out_stands
}


# Global functions ----
correct_IRRES2018_crowns <- function(data_trees, params_alloms, site) {
  
  data_trees_2018 <- data_trees %>% 
    dplyr::filter(year %in% 2015:2018) %>%
    
    dplyr::arrange(id_tree, desc(year)) %>% 
    
    dplyr::group_by(id_tree) %>%
    dplyr::mutate(
      
      dbh_cm_previous = dplyr::lead(dbh_cm),
      delta_dbh_cm = dbh_cm - dbh_cm_previous,
      delta_dbh_cm_positive = pmax(0, delta_dbh_cm),
      
      year_ref = ifelse(2015 %in% year, 2015, 2016),
      
      
      height_ref = ifelse(year_ref == 2015, 
                          h_m[year == 2015][1],
                          h_m[year == 2016][1]),
      height_was_measured = !is.na(height_ref),
      
      
      crownbase_ref = ifelse(year_ref == 2015, 
                             hbase_m[year == 2015][1],
                             hbase_m[year == 2016][1]),
      crownbase_was_measured = !is.na(crownbase_ref),
      
      
      radius_n_ref = ifelse(year_ref == 2015, 
                            rn_m[year == 2015][1],
                            rn_m[year == 2016][1]),
      radius_e_ref = ifelse(year_ref == 2015, 
                            re_m[year == 2015][1],
                            re_m[year == 2016][1]),
      radius_s_ref = ifelse(year_ref == 2015, 
                            rs_m[year == 2015][1],
                            rs_m[year == 2016][1]),
      radius_w_ref = ifelse(year_ref == 2015, 
                            rw_m[year == 2015][1],
                            rw_m[year == 2016][1]),
      
      radii_was_measured = 
        !is.na(radius_n_ref) & 
        !is.na(radius_e_ref) & 
        !is.na(radius_s_ref) & 
        !is.na(radius_w_ref)
      
    ) %>% 
    dplyr::ungroup() %>% 
    
    # Estimate delta height if height was measured in 2026
    dplyr::mutate(
      
      delta_height_m = ifelse(
        height_was_measured, 
        predict_height_increment(
          dbh_cm = dbh_cm_previous,
          delta_dbh_cm = delta_dbh_cm_positive,
          species = sp_allom,
          site = site,
          params_height = params_alloms$height
        ),
        0
      ),
      
      
      delta_radius_m = ifelse(
        radii_was_measured, 
        predict_crown_radius_increment(
          dbh_cm = dbh_cm_previous,
          delta_dbh_cm = delta_dbh_cm_positive,
          species = sp_allom,
          site = site,
          params_radius = params_alloms$cradius
        ),
        0
      )
      
    ) %>% 
    
    # Sum all increments from the current year up to 2026
    dplyr::group_by(id_tree) %>% 
    dplyr::mutate(
      height_increment_to_ref = rev(cumsum(coalesce(rev(delta_height_m), 0))),
      radius_increment_to_ref = rev(cumsum(coalesce(rev(delta_radius_m), 0)))
    ) %>% 
    dplyr::ungroup() %>% 
    
    dplyr::mutate(
      
      h_m = ifelse(
        
        !height_was_measured,
        
        # No measured 2015/2016 height : we will apply directly the static allometry
        NA,
        
        # Calculate next years dynamically
        height_ref + height_increment_to_ref
      ),
      
      
      rn_m = ifelse(
        
        !radii_was_measured,
        
        # No measured 2015/2016 radius : we will apply directly the static allometry
        NA,
        
        # Calculate next years
        radius_n_ref + radius_increment_to_ref
      ),
      
      re_m = ifelse(
        
        !radii_was_measured,
        
        # No measured 2015/2016 radius : we will apply directly the static allometry
        NA,
        
        # Calculate next years
        radius_e_ref + radius_increment_to_ref
      ),
      
      rs_m = ifelse(
        
        !radii_was_measured,
        
        # No measured 2015/2016 radius : we will apply directly the static allometry
        NA,
        
        # Calculate next years
        radius_s_ref + radius_increment_to_ref
      ),
      
      rw_m = ifelse(
        
        !radii_was_measured,
        
        # No measured 2015/2016 radius : we will apply directly the static allometry
        NA,
        
        # Calculate next years
        radius_w_ref + radius_increment_to_ref
      )
      
    ) %>% 
    dplyr::arrange(id_tree, desc(year)) %>% 
    
    dplyr::group_by(id_tree) %>%
    dplyr::mutate(
      
      height_m_previous = dplyr::lead(h_m),
      delta_height_m = h_m - height_m_previous,
      delta_height_m_positive = pmax(0, delta_height_m)
      
    ) %>% 
    dplyr::ungroup() %>% 
    
    # Estimate delta height if height was measured in 2026
    dplyr::mutate(
      
      delta_crownbase_m = ifelse(
        crownbase_was_measured, 
        predict_crownbase_increment(
          height_m = height_m_previous,
          delta_height_m = delta_height_m_positive,
          species = sp_allom,
          site = site,
          params_cbase = params_alloms$crownbase
        ),
        0
      )
      
    ) %>% 
    
    # Sum all increments from the current year up to 2026
    dplyr::group_by(id_tree) %>% 
    dplyr::mutate(crownbase_increment_to_ref = rev(cumsum(coalesce(rev(delta_crownbase_m), 0)))) %>% 
    dplyr::ungroup() %>% 
    
    dplyr::mutate(
      hbase_m = ifelse(
        
        !crownbase_was_measured,
        
        # No measured 2015/2016 crownbase : we will apply directly the static allometry
        NA,
        
        # Calculate next years
        crownbase_ref + crownbase_increment_to_ref
      )
    ) %>% 
    
    dplyr::filter(year == 2018)
  
  
  return(data_trees_2018)
}

estimate_missing_crowns <- function(data_trees, params_alloms, site) {
  
  data_trees %>% 
    
    dplyr::mutate(
      
      # Predict tree height
      height_measured = !is.na(h_m),
      
      h_m = ifelse(
        !height_measured,
        predict_height(dbh_cm, sp_allom, site, params_alloms$height),
        h_m
      ),
      
      
      # Predict tree crown radii
      radii_all_measured = !is.na(rn_m) & !is.na(re_m) & !is.na(rs_m) & !is.na(rw_m),
      
      rn_m = ifelse(
        is.na(rn_m),
        predict_crown_radius(dbh_cm, sp_allom, site, params_alloms$cradius),
        rn_m
      ),
      
      re_m = ifelse(
        is.na(re_m),
        predict_crown_radius(dbh_cm, sp_allom, site, params_alloms$cradius),
        re_m
      ),
      
      rs_m = ifelse(
        is.na(rs_m),
        predict_crown_radius(dbh_cm, sp_allom, site, params_alloms$cradius),
        rs_m
      ),
      
      rw_m = ifelse(
        is.na(rw_m),
        predict_crown_radius(dbh_cm, sp_allom, site, params_alloms$cradius),
        rw_m
      )
      
    ) %>% 
    
    dplyr::mutate(
      
      # Predict crown base height
      crownbase_measured = !is.na(hbase_m),
      
      hbase_m = ifelse(
        !crownbase_measured,
        predict_crownbase(h_m, sp_allom, site, params_alloms$crownbase),
        hbase_m
      )
      
    )
  
}

set_crown_shape <- function(data_trees, params_alloms) {
  
  data_trees %>% 
    
    dplyr::mutate(
      crown_type = ifelse(
        leaf_group == "broad", 
        "8E", # 8 eighth of ellipsoids for broadleaved
        "4P" # 4 fourth of paraboloids for needle-leaved
      ),
      
      hmax_measured = !is.na(hmax_m),
      
      hmax_m = case_when(
        !hmax_measured & crown_type == "8E" ~ hbase_m + 0.5 * (h_m - hbase_m),
        !hmax_measured & crown_type == "4P" ~ NA, # No need to specify it, it will automatically be setted to hbase_m in create_sl_stand()
        hmax_measured ~ hmax_m,
        TRUE ~ NA_real_
      )
      
    )
  
}

compute_volume_trees <- function(data_trees) {
  
  data_trees %>%
    
    # Specify if volume comes from full measurements (do not include hmax as a mandatory measurement)
    dplyr::mutate(
      volume_measured = height_measured & crownbase_measured & radii_all_measured
    ) %>% 
    
    # Estimate volume for either "4P" (fourth of paraboloids) or "8E" (eights of ellispoids)
    dplyr::mutate(
      
      # 8E: Upper and lower sections of ellipsoids
      c_top = dplyr::if_else(
        crown_type == "8E",
        h_m - hmax_m,
        NA_real_
      ),
      
      v_upper_n = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rn_m^2 * c_top / 8,
        NA_real_
      ),
      v_upper_s = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rs_m^2 * c_top / 8,
        NA_real_
      ),
      v_upper_w = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rw_m^2 * c_top / 8,
        NA_real_
      ),
      v_upper_e = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * re_m^2 * c_top / 8,
        NA_real_
      ),
      v_upper = v_upper_n + v_upper_s +
        v_upper_w + v_upper_e,
      
      c_bot = dplyr::if_else(
        crown_type == "8E",
        hmax_m - hbase_m,
        NA_real_
      ),
      
      v_lower_n = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rn_m^2 * c_bot / 8,
        NA_real_
      ),
      v_lower_s = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rs_m^2 * c_bot / 8,
        NA_real_
      ),
      v_lower_w = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * rw_m^2 * c_bot / 8,
        NA_real_
      ),
      v_lower_e = dplyr::if_else(
        crown_type == "8E",
        (4 / 3) * pi * re_m^2 * c_bot / 8,
        NA_real_
      ),
      v_lower = v_lower_n + v_lower_s +
        v_lower_w + v_lower_e,
      
      
      # 4P: Four quarters of paraboloids
      v_4p_n = dplyr::if_else(
        crown_type == "4P",
        pi * rn_m^2 * (h_m - hbase_m) / 8,
        NA_real_
      ),
      v_4p_s = dplyr::if_else(
        crown_type == "4P",
        pi * rs_m^2 * (h_m - hbase_m) / 8,
        NA_real_
      ),
      v_4p_w = dplyr::if_else(
        crown_type == "4P",
        pi * rw_m^2 * (h_m - hbase_m) / 8,
        NA_real_
      ),
      v_4p_e = dplyr::if_else(
        crown_type == "4P",
        pi * re_m^2 * (h_m - hbase_m) / 8,
        NA_real_
      ),
      
      # Total crown volume by crown type
      volume_m3 = dplyr::case_when(
        crown_type == "8E" ~ v_upper + v_lower,
        crown_type == "4P" ~ v_4p_n + v_4p_s +
          v_4p_w + v_4p_e,
        TRUE ~ NA_real_
      )
    )
}

compute_competition_trees <- function(data_trees, 
                                      local_ba_radius,
                                      plot_area, plot_xsize, plot_ysize) {
  
  # ── Plot weight (for per-hectare scaling) ──────────────────────────────────
  plot_weight <- 1 / plot_area
  local_weight <- 10000 / (pi * local_ba_radius^2)
  
  # ── Stand-level competition variables ─────────────────────────────────────
  data_trees <- data_trees %>% 
    dplyr::mutate(
      dg_cm         = sqrt(sum(dbh_cm^2) / n()),
      compet_status = dbh_cm / dg_cm,
      G_m2          = (pi * (dbh_cm / 100)^2 / 4),
      G_m2ha        = G_m2 * plot_weight,
      batot_m2ha    = sum(G_m2ha)
    ) %>%
    dplyr::arrange(dplyr::desc(dbh_cm)) %>%
    dplyr::mutate(bal_m2ha = ave(G_m2ha, FUN = cumsum) - G_m2ha) %>%
    dplyr::arrange(dplyr::desc(h_m)) %>%
    dplyr::mutate(bah_m2ha = ave(G_m2ha, FUN = cumsum) - G_m2ha)
  
  # ── Toroidal distance helper ───────────────────────────────────────────────
  # The shortest distance
  torus_dist <- function(xi, xj, size) {
    d <- xj - xi
    d - size * round(d / size)   # equivalent to modular wrap
  }
  
  # ── Local indices (computed row-by-row via sapply) ─────────────────────────
  n_trees <- nrow(data_trees)
  
  batot_local   <- numeric(n_trees)
  bal_local     <- numeric(n_trees)
  bah_local     <- numeric(n_trees)
  
  for (i in seq_len(n_trees)) {
    
    dx <- torus_dist(data_trees$x[i], data_trees$x, plot_xsize)
    dy <- torus_dist(data_trees$y[i], data_trees$y, plot_ysize)
    dist_to_i <- sqrt(dx^2 + dy^2)
    
    # Neighbours: within radius, excluding focal tree itself
    neighbours <- which(dist_to_i <= local_ba_radius & dist_to_i > 0)
    
    if (length(neighbours) == 0) {
      batot_local[i]  <- 0
      bal_local[i]    <- 0
      bah_local[i]    <- 0
    } else {
      G_neighbours     <- data_trees$G_m2[neighbours] * local_weight
      batot_local[i]   <- sum(G_neighbours)
      bal_local[i]     <- sum(G_neighbours[data_trees$dbh_cm[neighbours] > data_trees$dbh_cm[i]])
      bah_local[i]     <- sum(G_neighbours[data_trees$h_m[neighbours]   > data_trees$h_m[i]])
    }
  }
  
  # ── Attach local indices and return ───────────────────────────────────────
  data_trees %>%
    dplyr::mutate(
      batot_local_m2ha  = batot_local,
      bal_local_m2ha    = bal_local,
      bah_local_m2ha    = bah_local
    )
}



compute_light_variables <- function(sl_stand, data_rad_site) {
  
  tmp_out <- SamsaRaLight::run_sl(
    sl_stand, data_rad_site, 
    parallel_mode = TRUE,
    verbose = FALSE
  )
  
  sl_stand$sensors <- sl_stand$sensors %>% 
    dplyr::left_join(tmp_out$output$light$sensors[,c("id_sensor", "punobs")],
                     by = "id_sensor")
  
  sl_stand$trees <- sl_stand$trees %>% 
    dplyr::left_join(tmp_out$output$light$trees[,c("id_tree", "rci")],
                     by = "id_tree")
  
  return(sl_stand)
}


# FUNCTION ALLOMETRIES ----

predict_height <- function(dbh_cm, species, site, params_height) {
  
  data.frame(
    dbh_cm = dbh_cm,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_height$species, params_height$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_height, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      height_m = 1.3 + a1 * (1 - exp(-a2 * dbh_cm)) ^ a3
    ) %>% 
    dplyr::pull(height_m)
  
}

predict_height_increment <- function(dbh_cm, delta_dbh_cm, species, site, params_height) {
  
  data.frame(
    dbh_cm = dbh_cm,
    delta_dbh_cm = delta_dbh_cm,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_height$species, params_height$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_height, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      delta_height_m = a1 * a2 * a3 * exp( -a2 * dbh_cm) * (1 - exp(-a2 * dbh_cm)) ^ (a3 - 1) * delta_dbh_cm
    ) %>% 
    dplyr::pull(delta_height_m)
  
}

predict_crownbase <- function(height_m, species, site, params_cbase) {
  
  data.frame(
    height_m = height_m,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_cbase$species, params_cbase$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_cbase, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      crownbase_m = a * height_m ^ b
    ) %>% 
    dplyr::pull(crownbase_m)
  
}

predict_crownbase_increment <- function(height_m, delta_height_m, species, site, params_cbase) {
  
  data.frame(
    height_m = height_m,
    delta_height_m = delta_height_m,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_cbase$species, params_cbase$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_cbase, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      delta_crownbase_m = a * b * height_m ^ (b-1) * delta_height_m
    ) %>% 
    dplyr::pull(delta_crownbase_m)
  
}

predict_crown_radius <- function(dbh_cm, species, site, params_radius) {
  
  data.frame(
    dbh_cm = dbh_cm,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_radius$species, params_radius$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_radius, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      radius_m_mean = a * dbh_cm ^ b
    ) %>% 
    dplyr::pull(radius_m_mean)
  
}

predict_crown_radius_increment <- function(dbh_cm, delta_dbh_cm, species, site, params_radius) {
  
  data.frame(
    dbh_cm = dbh_cm,
    delta_dbh_cm = delta_dbh_cm,
    species = species,
    site = site
  ) %>% 
    
    # If no possible combination of species X site, set the population mean a1
    dplyr::mutate(
      site = dplyr::if_else(
        paste(species, site) %in% paste(params_radius$species, params_radius$site),
        site, "POPULATION"
      )
    ) %>% 
    
    dplyr::left_join(params_radius, by = c("species", "site")) %>% 
    
    dplyr::mutate(
      delta_radius_m_mean = a * b * dbh_cm ^ (b-1) * delta_dbh_cm
    ) %>% 
    dplyr::pull(delta_radius_m_mean)
  
}

