finalise_database <- function(init_db_raw,
                              data_plots) {

  init_db <- init_db_raw %>% 
    reorganise_database(data_plots) %>% 
    clean_database() %>% 
    add_species_info() %>% 
    plot_allometries(output_folder = "data/allometries_1_old", same_rad_color = TRUE) %>% 
    remove_allometries() %>% 
    plot_allometries(output_folder = "data/allometries_2_removed", same_rad_color = FALSE)
    

  return(init_db)
}



reorganise_database <- function(init_db_raw, data_plots) {
  
  # Get all unique site names across data types and origins
  site_names <- init_db_raw %>%
    purrr::map(~ purrr::map(.x, names)) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  # Reorganize the list by site
  purrr::map(site_names, function(site_name) {
    
    # Find the origin containing this site
    origins <- names(init_db_raw)[
      purrr::map_lgl(init_db_raw, function(origin_data) {
        any(purrr::map_lgl(origin_data, ~ site_name %in% names(.x)))
      })
    ]
    
    # Collect data types and origin names
    purrr::map(origins, function(origin_name) {
      
      origin_data <- init_db_raw[[origin_name]]
      
      # Add sites data
      site_data <- purrr::imap(origin_data, function(sites, data_type) {
        if (site_name %in% names(sites)) {
          sites[[site_name]]
        } else {
          NULL
        }
      }) %>%
        purrr::compact()
      
      # Add plot information as vector
      site_data$info <- data_plots %>% 
        dplyr::filter(site == site_name) %>% 
        dplyr::select(-site) %>% 
        as.vector()
      
      site_data
    }) %>%
      purrr::flatten()
  }) %>%
    rlang::set_names(site_names)
  
}


clean_database <- function(init_db) {
  
  # Clean for PACL between 0 and 1
  # Otherwise, wrong measurement ==> remove it and will be estimate with allometry
  init_db <- init_db %>%
    purrr::map(
      ~ .x %>%
        purrr::modify_in(
          "sensors",
          ~ .x %>%
            dplyr::filter(!is.na(pacl), pacl > 0, pacl < 1)
        )
    )
  
  
  # Hbase must be strictly lower than h
  # Hmax must be comprised between hbase and h
  # Otherwise, wrong measurement ==> remove it and will be estimate with allometry
  init_db <- init_db %>%
    purrr::map(
      ~ .x %>%
        purrr::modify_in(
          "trees",
          ~ .x %>%
            dplyr::mutate(
              invalid_height = !is.na(h_m) & !is.na(hbase_m) & h_m < hbase_m,
              invalid_hmax = !is.na(h_m) & !is.na(hbase_m) & !is.na(hmax_m) &
                (hmax_m > h_m | hmax_m < hbase_m),
              h_m = dplyr::if_else(invalid_height, NA_real_, h_m),
              hbase_m = dplyr::if_else(invalid_height, NA_real_, hbase_m),
              hmax_m = dplyr::if_else(invalid_hmax, NA_real_, hmax_m),
              rn_m = pmax(0.05, rn_m),
              re_m = pmax(0.05, re_m),
              rs_m = pmax(0.05, rs_m),
              rw_m = pmax(0.05, rw_m)
            ) %>%
            dplyr::select(-invalid_height, -invalid_hmax)
        )
    )
  
  # There is one tree with suspicious heights measurement
  # Seemed to be forgotten to convert in cm
  init_db$cloture_17$trees <- init_db$cloture_17$trees %>% 
    dplyr::mutate(
      h_m = ifelse(id_tree == "17_106", h_m * 100, h_m),
      hmax_m = ifelse(id_tree == "17_106", hmax_m * 100, hmax_m)
    )
  
  return(init_db)
}


add_species_info <- function(init_db) {
  
  # Observe the species in the database
  # species_by_plot <- init_db %>%
  #   purrr::map(
  #     ~ .x$trees %>%
  #       dplyr::distinct(species_name) %>%
  #       dplyr::arrange(species_name)
  #   ) %>%
  #   dplyr::bind_rows(.id = "site")
  # 
  # species_occ <- table(sort(
  #   init_db %>%
  #     purrr::map( ~ .x$trees) %>%
  #     dplyr::bind_rows() %>%
  #     dplyr::pull(species_name),
  #   decreasing = TRUE
  # ))
  
  # ABIES ALBA         ABIES GRANDIS   ACER PSEUDOPLATANUS               ACER sp       ALNUS GLUTINOSA       AUTRES FEUILLUS        BETULA PENDULA 
  # 213                     1                     1                    18                     2                     1                    11 
  # BETULA sp      BETULA VERRUCOSA           BROADLEAVES      CARPINUS BETULUS       CASTANEA SATIVA      CHENES INDIGENES       CORYLUS AVELANA 
  # 98                    38                    10                   279                     1                   601                    35 
  # CRATAEGUS sp       FAGUS SYLVATICA       ILEX AQUIFOLIUM          JUNIPERUS sp         LARIX DECIDUA              LARIX sp      MALUS SYLVESTRIS 
  # 1                  2582                     4                     1                   226                     1                     5 
  # PICEA ABIES      PICEA SITCHENSIS              PICEA sp      PINUS SYLVESTRIS       POPULUS TREMULA PSEUDOTSUGA MENZIESII       QUERCUS PETRAEA 
  # 5015                     3                    92                   904                    23                  3200                   990 
  # QUERCUS sp              SALIX sp      SORBUS AUCUPARIA             SORBUS sp     SORBUS TORMINALIS              TSUGA sp 
  # 2                     3                   371                    13                    38                    18 
  # 
  # Group species into specific calibration group 
  # that will be considered as modalities for the species random ffect
  sp_to_calib <- list(
    "Fagus_sylvatica" = "FAGUS SYLVATICA",
    "Quercus_sp" = c("QUERCUS sp", "QUERCUS PETRAEA", "CHENES INDIGENES"),
    "Carpinus_betulus" = "CARPINUS BETULUS",
    "Sorbus_sp" = c("SORBUS AUCUPARIA", "SORBUS sp", "SORBUS TORMINALIS"),
    "Populus_tremula" = "POPULUS TREMULA",
    "Betula_sp" = c("BETULA sp", "BETULA PENDULA", "BETULA VERRUCOSA"),
    "Acer_sp" = c("ACER PSEUDOPLATANUS", "ACER sp"),
    "Alnus_glutinosa" = "ALNUS GLUTINOSA",
    "Castanea_sativa" = "CASTANEA SATIVA", 
    "Corylus_avellana" = "CORYLUS AVELANA",
    "Crataegus_sp" = "CRATAEGUS sp", 
    "Ilex_aquifolium" = "ILEX AQUIFOLIUM",
    "Malus_sylvestris" = "MALUS SYLVESTRIS", 
    "Salix_sp" = "SALIX sp",
    "other_angio" = c("AUTRES FEUILLUS",  "BROADLEAVES"),
    
    "Picea_abies" = c("PICEA ABIES", "PICEA sp", "PICEA SITCHENSIS"),
    "Pseudotsuga_menziesii" = "PSEUDOTSUGA MENZIESII",
    "Pinus_sylvestris" = "PINUS SYLVESTRIS",
    "Abies_alba" = c("ABIES ALBA", "ABIES GRANDIS"),
    "Larix_decidua" = c("LARIX DECIDUA", "LARIX sp"),
    "Juniperus_sp" = "JUNIPERUS sp",
    "Tsuga_sp" = "TSUGA sp"
  )
  
  sp_to_calib_df <- tibble::tibble(
    sp_in_inv = sp_to_calib,
    sp_calib = names(sp_to_calib)
  ) %>%
    tidyr::unnest(sp_in_inv) %>%
    
    dplyr::mutate(
      
      # Set the reference species for fitting allometries to reconstruct the 3D stands
      sp_allom = dplyr::case_when(
        
        # Angiosperms with sufficient occurrences
        sp_calib %in% c(
          "Fagus_sylvatica",
          "Quercus_sp",
          "Carpinus_betulus",
          "Sorbus_sp",
          "Betula_sp"
        ) ~ sp_calib,
        
        # All remaining angiosperms
        sp_calib %in% c(
          "Acer_sp",
          "Alnus_glutinosa",
          "Castanea_sativa",
          "Corylus_avellana",
          "Crataegus_sp",
          "Ilex_aquifolium",
          "Malus_sylvestris",
          "Populus_tremula",
          "Salix_sp",
          "other_angio"
        ) ~ "other_angio",
        
        # Gymnosperms: assign rare species to proxy allometries
        sp_calib %in% c("Abies_alba", "Tsuga_sp") ~ "Abies_alba",
        
        sp_calib %in% c(
          "Pinus_sylvestris",
          "Juniperus_sp"
        ) ~ "Pinus_sylvestris",
        
        # Other gymnosperms
        sp_calib %in% c(
          "Picea_abies",
          "Pseudotsuga_menziesii",
          "Larix_decidua"
        ) ~ sp_calib,
        
        TRUE ~ NA_character_
      ),
      
      # Species that contribute to fitting
      use_for_allom = !sp_calib %in% c(
        "Tsuga_sp", # Use Abies alba allometry for rpediction without using Tsuga observaitons for fitting
        "Juniperus_sp" # Use Pinus sylvestris allometry for rpediction without using Tsuga observaitons for fitting
      ),
      
      # Set the species-specific shade tolerance value
      shadetol = dplyr::case_match(
        sp_calib,
        
        "Fagus_sylvatica" ~ 4.56,
        "Quercus_sp" ~ mean(c(2.73, 2.45)), # Quercus petraea, Quercus robur
        "Carpinus_betulus" ~ 3.97,
        "Sorbus_sp" ~ 2.73, # Sorbus aucuparia (it is greatly dominant in the dataset, 371/422: 88%)
        "Populus_tremula" ~ 2.22,
        "Betula_sp" ~ mean(c(2.03, 1.85)), # Betula pendula, Betula pubescens
        "Acer_sp" ~ mean(c(4.2, 3.73)), # Acer platanoides, Acer pseudoplatanus
        "Alnus_glutinosa" ~ 2.71,
        "Castanea_sativa" ~ 3.15, 
        "Corylus_avellana" ~ 3.53,
        "Crataegus_sp" ~ mean(c(2.45, 1.93)), # Crataegus laevigata agg., Crataegus monogyna agg.
        "Ilex_aquifolium" ~ 3.86,
        "Malus_sylvestris" ~ 2.32, 
        "Salix_sp" ~ 1.99, # Salix alba by default (3 individuals in cloture 18 only)
        "other_angio" ~ 2.45, # Mean of all angio in the Niinemets and Valladares database
        
        "Picea_abies" ~ 4.45,
        "Pseudotsuga_menziesii" ~ 2.78,
        "Pinus_sylvestris" ~ 1.67,
        "Abies_alba" ~ 4.6,
        "Larix_decidua" ~ 1.46,
        "Juniperus_sp" ~ 1.71, # Juniperus communis by default (1 individual in Ventoux_C7 only)
        "Tsuga_sp" ~ 4.96
      ),
      
      # Set the species leaf type
      leaf_group = dplyr::case_match(
        sp_calib,
        c(
          "Fagus_sylvatica", "Quercus_sp", "Carpinus_betulus",
          "Sorbus_sp", "Populus_tremula", "Betula_sp", 
          "Acer_sp", "Alnus_glutinosa", "Castanea_sativa", 
          "Corylus_avellana", "Crataegus_sp", "Ilex_aquifolium",
          "Malus_sylvestris", "Salix_sp", "other_angio"
        ) ~ "broad",
        c(
          "Picea_abies", "Pseudotsuga_menziesii",
          "Pinus_sylvestris", "Abies_alba",
          "Larix_decidua", "Juniperus_sp", "Tsuga_sp"
        ) ~ "needle"
      )
    )
  
  
  # Add to inventories
  init_db <- init_db %>%
    purrr::map(
      ~ .x %>% purrr::modify_in(
        "trees",
        ~ .x %>%
          dplyr::left_join(
            sp_to_calib_df,
            by = c("species_name" = "sp_in_inv")
          ) %>% 
          dplyr::rename(sp_in_inv = species_name)
      )
    )
  
  return(init_db)
}


plot_allometries <- function(init_db, output_folder, same_rad_color = TRUE) {
  
  for (i in seq_along(init_db)) {
    
    site_name <- names(init_db)[i]
    origin <- init_db[[i]]$info$origin
    
    plot_folder <- file.path(output_folder, origin)
    dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
    
    init_db[[i]]$trees <- init_db[[i]]$trees %>%
      mutate(
        crown_ratio = if_else(
          is.finite(h_m) & h_m > 0 & is.finite(hbase_m),
          1 - hbase_m / h_m,
          NA_real_
        ),
        same_rad = case_when(
          is.na(rn_m) | is.na(re_m) | is.na(rs_m) | is.na(rw_m) ~ NA,
          TRUE ~ abs(rn_m - re_m) < 1e-6 &
            abs(rn_m - rs_m) < 1e-6 &
            abs(rn_m - rw_m) < 1e-6
        )
      )
    
    plot_specs <- list(
      H_vs_DBH = list(
        x = "dbh_cm", y = "h_m",
        xlab = "DBH (cm)", ylab = "Height (m)",
        title = "Height ~ DBH"
      ),
      RN_vs_DBH = list(
        x = "dbh_cm", y = "rn_m",
        xlab = "DBH (cm)", ylab = "North crown radius (m)",
        title = "RN ~ DBH"
      ),
      CBH_vs_H = list(
        x = "h_m", y = "hbase_m",
        xlab = "Height (m)", ylab = "Crown base height (m)",
        title = "CBH ~ H"
      ),
      HMAX_vs_H = list(
        x = "h_m", y = "hmax_m",
        xlab = "Height (m)", ylab = "Crown maximum height (m)",
        title = "HMAX ~ H"
      ),
      crown_ratio_vs_DBH = list(
        x = "dbh_cm", y = "crown_ratio",
        xlab = "DBH (cm)", ylab = "Crown ratio",
        title = "Crown ratio ~ DBH"
      )
    )
    
    for (plot_name in names(plot_specs)) {
      
      spec <- plot_specs[[plot_name]]
      
      # Exclude missing/non-finite x and y values to avoid geom_point warnings
      plot_data <- init_db[[i]]$trees %>%
        filter(
          is.finite(.data[[spec$x]]),
          is.finite(.data[[spec$y]])
        )
      
      p <- ggplot(
        plot_data,
        aes(
          x = .data[[spec$x]],
          y = .data[[spec$y]],
          color = if (same_rad_color) factor(same_rad) else factor(sp_in_inv)
        )
      ) +
        geom_point(alpha = 0.7, na.rm = TRUE) +
        labs(
          x = spec$xlab,
          y = spec$ylab,
          color = if (same_rad_color) "Same radius" else "Species",
          title = paste(spec$title, "-", origin, "-", site_name)
        ) +
        theme_bw()
      
      ggsave(
        filename = file.path(
          plot_folder,
          paste0(plot_name, "_", site_name, ".png")
        ),
        plot = p,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
    
  }
  
  return(init_db)
}


remove_allometries <- function(init_db) {
  
  # For IRRES, there is no allometries, all values have been measured
  # For INRAe plots: Ventoux, Lorris and Frontdefrey, no allometries, all values have been measured
  
  # For Cloture: height and crown base heights values all comes from measurements
  # Except for few trees in coming from missing measurements
  # To detect them, we associate a value with more than 1 decimal to an allometric relationships
  # As the measured values are more than 1 decimal of precision (standrad in cloture for measured values)
  # Whereas allometric values come from equations, thus with several decimals
  # And remove radisu measurements when the four directions are equals ==> allometries
  init_db <- init_db %>%
    purrr::imap(~ {
      if (grepl("cloture", .y)) {
        .x$trees <- .x$trees %>%
          dplyr::mutate(
            
            remove_h = abs(h_m * 10 - round(h_m * 10)) > 1e-8,
            remove_hbase = abs(hbase_m * 10 - round(hbase_m * 10)) > 1e-8,
            remove_heights = remove_h & remove_hbase, 
            # Both are in large decimals = high certainty that it is allometric, because height and cbh are always taken together
              
            h_m = ifelse(remove_heights, NA, h_m),
            hbase_m = ifelse(remove_heights, NA, hbase_m),

            rn_m = ifelse(same_rad, NA, rn_m),
            re_m = ifelse(same_rad, NA, re_m),
            rs_m = ifelse(same_rad, NA, rs_m),
            rw_m = ifelse(same_rad, NA, rw_m)
          ) %>% 
          dplyr::select(-remove_h, -remove_hbase, -remove_heights)
      }
      .x
    })
  
  
  # For Baileux plots
  init_db <- init_db %>%
    purrr::imap(~ {
      if (grepl("baileux", .y)) {
        .x$trees <- .x$trees %>%
          dplyr::mutate(
            
            # Allometric relationships for height if (DBH > 20cm and same radius and oak stand) that seems to be allometric relationships
            remove_h = same_rad & .y == "baileux_oak" & dbh_cm > 20,
            
            # Allometric relationships for crown ratio if same radius and (DBH < 20cm and beech stand) that seems to be allometric relationships
            # And we remove all allometries for CBH in oak stand as it is too difficult to discriminate alloemtries and they are in large majority
            remove_hbase = same_rad | (.y == "baileux_beech" & dbh_cm < 20) | (.y == "baileux_oak"),

            h_m = ifelse(remove_h, NA, h_m),
            hbase_m = ifelse(remove_hbase, NA, hbase_m),
            hmax_m = ifelse(remove_hbase, NA, hmax_m),
            
            # And remove radius measurements when the four directions are equals ==> allometries
            rn_m = ifelse(same_rad, NA, rn_m),
            re_m = ifelse(same_rad, NA, re_m),
            rs_m = ifelse(same_rad, NA, rs_m),
            rw_m = ifelse(same_rad, NA, rw_m)
          ) %>% 
          dplyr::select(-remove_h, -remove_hbase)
      }
      .x
    })
  
  return(init_db)
}


