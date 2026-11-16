format_database <- function(data_irres_fp, 
                            data_cloture_fp, 
                            data_baileux_folderpath,
                            data_fontfreyde_folderpath,
                            data_lorris_folderpath,
                            data_ventoux_folderpath,
                            dico_species_fp) {
  
  data_irres2018 <- format_data_irres(data_irres_fp, dico_species_fp, years = c(2015, 2016, 2018), prefix = "irres2018")
  data_irres2026 <- format_data_irres(data_irres_fp, dico_species_fp, years = 2026, prefix = "irres2026")
  data_cloture <- format_data_cloture(data_cloture_fp, dico_species_fp)
  data_baileux <- format_data_baileux(data_baileux_folderpath, dico_species_fp)
  data_fontfreyde <- format_data_fontfreyde(data_fontfreyde_folderpath, dico_species_fp)
  data_lorris <- format_data_lorris(data_lorris_folderpath, dico_species_fp)
  data_ventoux <- format_data_ventoux(data_ventoux_folderpath, dico_species_fp)
  
  return(
    list(
      "irres2018" = data_irres2018,
      "irres2026" = data_irres2026,
      "cloture" = data_cloture,
      "baileux" = data_baileux,
      "fontfreyde" = data_fontfreyde,
      "lorris" = data_lorris,
      "ventoux" = data_ventoux
    )
  )
}


# Formatting functions ----
format_data_fontfreyde <- function(data_fontfreyde_folderpath, dico_species_fp) {
  
  message("----- Formatting fontfreyde database...")
  
  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "PB")

  # Create tree dataset
  id_plots <- 1:5
  data_trees <- setNames(vector("list", length(id_plots)), paste0("fontfreyde_", id_plots))
    
  for(i in id_plots){
    
    tmp_id_plot <- paste0("fontfreyde_", i) 
    
    tmp_stand <- readxl::read_excel(
      file.path(
        data_fontfreyde_folderpath,
        paste0("Fontfreyde_P", i, "_Houppier_ApEclaircie.xls")
      )
    ) %>%
      dplyr::mutate(
        id_plot = tmp_id_plot,
        year = 2007,
        hmax_m = NA,
        X = as.numeric(X) / 100,
        Y = as.numeric(Y) / 100
      ) %>%
      dplyr::rename(
        id_tree_str = NArbre,
        species_name = Essence,
        h_m = H,
        hbase_m = `Base Houp`,
        dbh_cm = DBH,
        x = X,
        y = Y
      ) %>%
      dplyr::mutate(
        id_tree = dplyr::row_number()
      ) %>% 
      
      # Remove NA radius that corresponds to multistems. 
      # At each multistem, one stem has been measures considering the whole crown
      tidyr::drop_na(r_houp1) %>% 
      
      # Harmonize species name
      dplyr::left_join(dico_sp[,c("PB_species", "IRRES_code")], by = c("species_name" = "PB_species"))     
    
    
      # Crown measurements
    tmp_radius <- tmp_stand %>%
      
      # Get radius and azimuth columns
      dplyr::select(
        id_tree, id_tree_str,
        starts_with("r_houp"),
        starts_with("a_houp")
      ) %>%
      
      # Convert r_houp1 / a_houp1 etc. into rows
      tidyr::pivot_longer(
        cols = -c(id_tree, id_tree_str),
        names_to = c(".value", "id_radius"),
        names_pattern = "([ra])_houp(.*)"
      ) %>%
      dplyr::rename(
        radius_cm = r,
        azimuth = a
      ) %>%
      tidyr::drop_na(radius_cm, azimuth) %>% 
      dplyr::mutate(
        radius_m = radius_cm / 100
      ) %>% 
      
      # Assign the closest unique cardinal
      assign_closest_cardinal() %>%
      
      tidyr::pivot_wider(
        id_cols = c(id_tree, id_tree_str),
        names_from = cardinal,
        values_from = radius_m,
        names_glue = "r{cardinal}_m"
      )
      
    
    data_trees[[tmp_id_plot]] <- tmp_stand %>% 
      dplyr::left_join(tmp_radius, by = "id_tree") %>% 
      dplyr::select(
        id_tree, year,
        species_name,
        x, y, 
        dbh_cm, h_m, hbase_m, hmax_m, 
        rn_m, re_m, rs_m, rw_m
      ) %>% 
      dplyr::arrange(id_tree, year)
  }
  
  
  # Create data sensors
  data_sensors <- readxl::read_excel(
    file.path(data_fontfreyde_folderpath, "Fontfreyde2007_Transmis.xls")
  ) %>%
    dplyr::mutate(
      id_plot = paste0("fontfreyde_", Placette),
      id_sensor = row_number(),
      year = 2007,
      x = as.numeric(X) / 100,
      y = as.numeric(Y) / 100,
      h_m = 1.3,
      pacl = as.numeric(Transmis_cor)
    ) %>% 
    dplyr::select(
      id_plot, id_sensor, year,
      x, y, h_m, pacl
    ) %>%

    # Group into a list of trees by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot))  
  
  
  # Return data IRRES
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data.frame() # No information here, we will use automatic detection
    )
  )
  
}


format_data_lorris <- function(data_lorris_folderpath, dico_species_fp) {
  
  message("----- Formatting lorris database...")
  
  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "PB")
  
  # Create tree dataset
  id_plots <- c(16, 38, 255)
  data_trees <- setNames(vector("list", length(id_plots)), paste0("lorris_", id_plots))

  for(i in id_plots){
    
    tmp_id_plot <- paste0("lorris_", i) 
    
    tmp_stand <- readxl::read_excel(
      file.path(
        data_lorris_folderpath,
        paste0("Lorris_Parcelle", i, "_Final.xls")
      ),
      sheet = 2
    ) %>%
      dplyr::mutate(
        id_plot = tmp_id_plot,
        year = 2007,
        hmax_m = NA
      ) %>%
      dplyr::rename(
        id_tree_str = id_arbre,
        species_name = essence,
        h_m = H,
        hbase_m = `Base Houp`,
        dbh_cm = DBH,
        x = X_dispo,
        y = Y_dispo
      ) %>%
      dplyr::mutate(
        id_tree = dplyr::row_number()
      ) %>% 
      
      # Remove NA radius that corresponds to multistems or one DEAD tree in Lorris1
      # At each multistem, one stem has been measures considering the whole crown
      tidyr::drop_na(r_houp1) %>% 
      
      # Harmonize species name
      dplyr::left_join(dico_sp[,c("PB_species", "IRRES_code")], by = c("species_name" = "PB_species"))     
    
    
    # Crown measurements
    tmp_radius <- tmp_stand %>%
      
      # Get radius and azimuth columns
      dplyr::select(
        id_tree, id_tree_str,
        starts_with("r_houp"),
        starts_with("a_houp")
      ) %>%
      
      # Convert r_houp1 / a_houp1 etc. into rows
      tidyr::pivot_longer(
        cols = -c(id_tree, id_tree_str),
        names_to = c(".value", "id_radius"),
        names_pattern = "([ra])_houp(.*)"
      ) %>%
      dplyr::rename(
        radius_cm = r,
        azimuth = a
      ) %>%
      tidyr::drop_na(radius_cm, azimuth) %>% 
      dplyr::mutate(
        radius_m = radius_cm / 100
      ) %>% 
      
      # Assign the closest unique cardinal
      assign_closest_cardinal() %>%
      
      tidyr::pivot_wider(
        id_cols = c(id_tree, id_tree_str),
        names_from = cardinal,
        values_from = radius_m,
        names_glue = "r{cardinal}_m"
      )
    
    
    data_trees[[tmp_id_plot]] <- tmp_stand %>% 
      dplyr::left_join(tmp_radius, by = "id_tree") %>% 
      dplyr::select(
        id_tree, year,
        species_name,
        x, y, 
        dbh_cm, h_m, hbase_m, hmax_m, 
        rn_m, re_m, rs_m, rw_m
      ) %>% 
      dplyr::arrange(id_tree, year)
    
  }
  
  
  # Create data sensors
  data_sensors <- vroom::vroom(
    file.path(
      data_lorris_folderpath,
      "data_pacl_lorris.csv"
    ),
    show_col_types = F
  ) %>%
    dplyr::mutate(
      id_plot = paste0("lorris_", id_plot),
      id_sensor = row_number(),
      year = 2007
    ) %>% 
    dplyr::select(
      id_plot, id_sensor, year,
      x, y, h_m = z, pacl = PACLtotal
    ) %>%
    
    # Group into a list of sensors by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot))
  
  
  # Return data IRRES
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data.frame() # No information here, we will use automatic detection
    )
  )
  
}


format_data_ventoux <- function(data_ventoux_folderpath, dico_species_fp) {
  
  message("----- Formatting ventoux database...")
  
  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "PB")
  
  # Create tree dataset
  id_plots <- c("34", "C7")
  data_trees <- setNames(vector("list", length(id_plots)), paste0("ventoux_", id_plots))
  
  for(i in id_plots){
    
    tmp_id_plot <- paste0("ventoux_", i) 
    
    tmp_stand <- readxl::read_excel(
      file.path(
        data_ventoux_folderpath,
        paste0("Ventoux_placette", i, "_Houppier.xls")
      ),
      sheet = 2
    ) %>%
      dplyr::mutate(
        id_plot = tmp_id_plot,
        year = 2008,
        hmax_m = NA,
        hauteur = hauteur / 100,
        h_houppier = h_houppier / 100,
        X = as.numeric(X),
        Y = as.numeric(Y)
      ) %>%
      dplyr::rename(
        id_tree_str = id_arbre,
        species_name = Essence,
        h_m = hauteur,
        hbase_m = h_houppier,
        dbh_cm = DBH,
        x = X,
        y = Y
      ) %>%
      dplyr::mutate(
        id_tree = dplyr::row_number()
      ) %>% 
      
      # Remove NA heights that corresponds to multistems. 
      # At each multistem, one stem has been measures considering the whole crown
      tidyr::drop_na(houpR1) %>% 
      
      # Harmonize species name
      dplyr::left_join(dico_sp[,c("PB_species", "IRRES_code")], by = c("species_name" = "PB_species"))     
    
    
    # Trees' position are in relative cm only in Ventoux_C7
    if (i == "C7") {
      tmp_stand$x <- tmp_stand$x / 100
      tmp_stand$y <- tmp_stand$y / 100
    }
    
    
    # Crown measurements
    tmp_radius <- tmp_stand %>%
      
      # Get radius and azimuth columns
      dplyr::select(
        id_tree, id_tree_str,
        starts_with("houpR"),
        starts_with("houpAz")
      ) %>% 
      
      # Convert r_houp1 / a_houp1 etc. into rows
      tidyr::pivot_longer(
        cols = -c(id_tree, id_tree_str),
        names_to = c(".value", "id_radius"),
        names_pattern = "houp(R|Az)(.*)"
      ) %>%
      
      dplyr::rename(
        radius_cm = R,
        azimuth = Az
      ) %>%
      tidyr::drop_na(radius_cm, azimuth) %>% 
      dplyr::mutate(
        radius_m = radius_cm / 100
      ) %>% 
      
      # Assign the closest unique cardinal
      assign_closest_cardinal() %>%
      
      tidyr::pivot_wider(
        id_cols = c(id_tree, id_tree_str),
        names_from = cardinal,
        values_from = radius_m,
        names_glue = "r{cardinal}_m"
      )
    
    
    data_trees[[tmp_id_plot]] <- tmp_stand %>% 
      dplyr::left_join(tmp_radius, by = "id_tree") %>% 
      dplyr::select(
        id_tree, year,
        species_name,
        x, y, 
        dbh_cm, h_m, hbase_m, hmax_m, 
        rn_m, re_m, rs_m, rw_m
      ) %>% 
      dplyr::arrange(id_tree, year)
    
  }
  
  
  # Create data sensors
  data_sensors <- vroom::vroom(
    file.path(
      data_ventoux_folderpath,
      "data_pacl_ventoux.csv"
    ),
    show_col_types = F
  ) %>%
    dplyr::mutate(
      id_plot = paste0("ventoux_", id_plot),
      id_sensor = row_number(),
      year = 2008
    ) %>% 
    dplyr::select(
      id_plot, id_sensor, year,
      x, y, h_m = z, pacl = PACLtotal
    ) %>%

    # Group into a list of sensors by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot))
  
  
  # Return data IRRES
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data.frame() # No information here, we will use automatic detection
    )
  )
  
}


format_data_baileux <- function(data_baileux_folderpath, dico_species_fp) {
  
  message("----- Formatting baileux database...")
  
  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "Baileux")
  
  # Start end end lines in the .txt file for the different sections
  df <- data.frame(
    type = c("beech", "mixed", "oak"),
    start_tree = c(45, 45, 42),
    end_tree = c(512, 593, 437),
    start_sensor = c(18, 18, 18),
    end_sensor = c(28, 26, 22),
    start_polygon = c(537, 617, 461),
    end_polygon = c(540, 620, 464)
  )
  
  # For each site
  data_trees <- setNames(vector("list", nrow(df)), df$type)
  data_sensors <- setNames(vector("list", nrow(df)), df$type)
  data_polygons <- setNames(vector("list", nrow(df)), df$type)
  
  for(i in seq_len(nrow(df))){
    
    tmp_fp <- file.path(
      data_baileux_folderpath, 
      paste0("Baileux_", df$type[i],"_2012_extended2003.xlsx")
    )
    
    # tree data frame
    data_trees[[df$type[i]]] <- suppressMessages(readxl::read_xlsx(tmp_fp)) %>% 
      
      # Transform and clean tree dataset from text file
      dplyr::slice(df$start_tree[i]:df$end_tree[i]) %>%
      janitor::row_to_names(row_number = 1) %>%
      dplyr::select(
        id = `#id`,
        speciesId = speciesName,
        h_m = `htot (m)`,
        hmax_m = `hlce (m)`,
        hbase_m = `hcb (m)`,
        rn_m = `rnorth (m)`,
        re_m = `reast (m)`,
        rs_m = `rsouth (m)`,
        rw_m = `rwest (m)`,
        gbh_cm = `gbh (cm)`,
        x, y
      ) %>%
      replace((.) == "NA", NA) %>%
      dplyr::mutate(
        id_tree = as.numeric(id),
        year = 2012,
        speciesId = as.numeric(speciesId),
        dbh_cm = as.numeric(gbh_cm)/pi,
        rn_m = as.numeric(rn_m),
        re_m = as.numeric(re_m),
        rs_m = as.numeric(rs_m),
        rw_m = as.numeric(rw_m),
        x = as.numeric(x),
        y = as.numeric(y),
        h_m = as.numeric(h_m),
        hbase_m = as.numeric(hbase_m),
        hmax_m = as.numeric(hmax_m),
        same_rad = pmin(rn_m, re_m, rs_m, rw_m) == pmax(rn_m, re_m, rs_m, rw_m)
      ) %>% 
      dplyr::mutate_if(is.numeric, round, 4) %>% 
      
      # Add the species
      dplyr::left_join(dico_sp[,c("niceName", "speciesId")], by = c("speciesId")) %>% 
      dplyr::mutate(niceName = toupper(niceName)) %>% 
      
      # Format variables
      dplyr::select(
        id_tree, year,
        species_name = niceName,
        x, y, 
        dbh_cm, h_m, hbase_m, hmax_m, 
        rn_m, re_m, rs_m, rw_m
      ) %>% 
      dplyr::arrange(id_tree, year)
    
    
    # sensors data.frame
    data_sensors[[df$type[i]]] <- suppressMessages(readxl::read_xlsx(tmp_fp))[1:7] %>%
      dplyr::slice(df$start_sensor[i]:df$end_sensor[i]) %>%
      janitor::row_to_names(row_number = 1) %>%
      dplyr::rename(
        ID = `# ID`,
        z = h
      ) %>%
      dplyr::mutate(
        id_sensor = row_number(),
        year = 2012,
        x = as.numeric(x),
        y = as.numeric(y),
        h_m = as.numeric(z),
        pacl = as.numeric(PACLTotal)
      ) %>%
      dplyr::select(
        id_sensor, year, 
        x, y, h_m, pacl
      )
    
    
    # plot extents data.frame
    data_polygons[[df$type[i]]] <- suppressMessages(readxl::read_xlsx(tmp_fp))[2:3] %>%
      dplyr::slice(df$start_polygon[i]:df$end_polygon[i]) %>%
      dplyr::rename(x = `...2`, y = `...3`) %>% 
      dplyr::mutate(
        x = as.numeric(x),
        y = as.numeric(y)
      ) %>% 
      dplyr::select(x, y)
  }
  
  names(data_trees) = paste0("baileux_", names(data_trees))
  names(data_sensors) = paste0("baileux_", names(data_sensors))
  names(data_polygons) = paste0("baileux_", names(data_polygons))
  
  # Return data Cloture
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data_polygons
    )
  )
}


format_data_cloture <- function(data_cloture_fp, dico_species_fp) {
  
  message("----- Formatting cloture database...")
  
  load(data_cloture_fp)

  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "cloture")
  
  # Format tree inventory
  data_trees <- trees %>%
    
    # Keep only the first observed year (year of crown and dbh measurement)
    dplyr::group_by(plot, field_plot, field_number) %>%
    dplyr::arrange(year, by_group = TRUE) %>% 
    dplyr::mutate(
      n_measures = n(),
      id_measure = row_number()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(id_measure == 1) %>% 
    
    # Create ids
    dplyr::mutate(
      id_tree = row_number(),
      id_plot = paste0("cloture_", plot),
      hmax_m = NA
    ) %>% 
    
    # Harmonize species name
    dplyr::left_join(dico_sp[,c("species", "latin_name")], by = "species") %>% 
    
    # Format variable names
    dplyr::select(
      id_plot, plot, subplot = field_plot, id_tree, year,
      species_name = latin_name, 
      x, y, dbh_cm = dbh, 
      h_m = height, hbase_m = crownBaseHeight, hmax_m, 
      rn_m = RN, re_m = RE, rs_m = RS, rw_m = RW
    ) %>% 
    dplyr::arrange(id_plot, id_tree, year) %>% 
    
    # Group into a list of trees by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot)) 
  
  
  # Get data radiation
  data_sensors <- data_trees %>% 
    dplyr::bind_rows(.id = "id_plot") %>% 
    dplyr::select(id_plot, plot, year) %>% 
    dplyr::distinct() %>% 
    dplyr::left_join(sensors, by = c("plot", "year")) %>% 
    
    dplyr::select(
      id_plot, id_sensor = subplot_id, 
      year,
      x, y, h_m = h, 
      pacl = pacl_total
    ) %>% 
    
    # Group into a list of trees by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot)) 
  
  
  
  # Get plot extents
  # Initial inventory ==> rem == R2
  # But some plots had an extended inventory, with three different area depending on trees DBH
  # E2 < E3 < E4
  
  # > table(plot_extents$plot, plot_extents$rem)
  # 
  # E2  E3  E4  R2
  # 1   42  57 101  55
  # 2    0   0   0  55
  # 3   38  55 101  53
  # 4   40  58 101  56
  # 5    0   0   0  54
  # 8    0   0   0  55
  # 9   40  56 103  56
  # 10  40  57 100  56
  # 11  40  57 101  55
  # 12   0   0   0  55
  # 14   0   0   0  55
  # 15   0   0   0  56
  # 17  36  54  96  51
  # 18   0   0   0 114
  # 19   0   0   0  98
  # 20   0   0   0  94
  # 22  38  55  99  54
  # 23  38  55  98  52
  # 24  37  54  98  52
  # 25   0   0   0 149
  # 26   0   0   0  56
  # 28   0   0   0  56
  # 29   0   0   0  55
  
  # for (i in unique(plot_extents$plot)) {
  #   ggsave(
  # 
  #     file.path("data/cloture", paste0("cloture_", i, ".png")),
  # 
  #     ggplot() +
  #       # Plot polygons, colored by remark
  #       geom_polygon(
  #         data = plot_extents %>% dplyr::filter(plot == i),
  #         aes(x = x, y = y, fill = rem),
  #         color = "black",
  #         alpha = 0.5
  #       ) +
  #       # Tree locations
  #       geom_point(
  #         data = data_trees[[paste0("cloture_", i)]],
  #         aes(x = x, y = y),
  #         color = "black",
  #         size = 1,
  #         alpha = 0.8
  #       ) +
  #       coord_equal() +
  #       labs(
  #         x = "X",
  #         y = "Y",
  #         fill = "Remark",
  #         title = paste("Plot", i)
  #       ) +
  #       theme_minimal()
  #   )
  # }

  
  data_polygons <- plot_extents %>% 
    dplyr::group_by(plot) %>% 
    dplyr::filter(if("E4" %in% rem) rem == "E4" else rem == "R2") %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(id_plot = paste0("cloture_", plot)) %>% 
    dplyr::select(-plot) %>% 
    
    # Group into a list of trees by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot)) 
  
  
  # Remove trees outside the polygon
  data_trees <- purrr::map2(
    data_trees,
    data_polygons[names(data_trees)],
    function(trees, polygon) {
      
      # Polygon vertices
      poly_x <- polygon$x
      poly_y <- polygon$y
      
      # Remove missing polygon coordinates
      valid <- is.finite(poly_x) & is.finite(poly_y)
      poly_x <- poly_x[valid]
      poly_y <- poly_y[valid]
      
      # Close polygon if necessary
      if (length(poly_x) > 0 &&
          (poly_x[1] != poly_x[length(poly_x)] ||
           poly_y[1] != poly_y[length(poly_y)])) {
        poly_x <- c(poly_x, poly_x[1])
        poly_y <- c(poly_y, poly_y[1])
      }
      
      # Check which trees are inside or on the boundary
      inside <- sp::point.in.polygon(
        point.x = trees$x,
        point.y = trees$y,
        pol.x = poly_x,
        pol.y = poly_y
      )
      
      # Values: 0 = outside, 1 = inside, 2 or 3 = boundary
      dplyr::filter(trees, inside > 0)
    }
  )
  
  
  # Return data Cloture
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data_polygons
    )
  )
}

  
format_data_irres <- function(data_irres_fp, dico_species_fp, years, prefix) {
  
  message("----- Formatting ", prefix, " database...")
  
  load(data_irres_fp)
  
  # Get species dico
  dico_sp <- readxl::read_excel(dico_species_fp, sheet = "IRRES")
  
  # Format tree inventory
  # Here we start from the raw complete IRRES database from 2015 to 2026
  data_trees <- db %>% 
    get_irres_trees() %>% 
    check_timeseries_continuity() %>%
    clean_discontinuities() %>% 
    check_timeseries_continuity() %>%
    estimate_dbh_2015() %>% 
    keep_living_trees() %>% 
    estimate_missing_dbh() %>% 
    keep_adult_trees() %>% 
    add_crown_measurements(db) %>% 
    keep_years(years) %>% 
    
    # Harmonize species name
    dplyr::left_join(dico_sp[,c("Code", "Latin_name")], by = c("speciesId" = "Code")) %>% 
    
    # Format variable names
    dplyr::mutate(
      id_plot = paste0(prefix, "_", id_plot),
      hmax_m = NA
    ) %>% 
    dplyr::select(
      id_plot, id_tree, year,
      species_name = Latin_name,
      x, y, 
      dbh_cm, h_m, hbase_m, hmax_m, 
      rn_m, re_m, rs_m, rw_m
    ) %>% 
    dplyr::arrange(id_plot, id_tree, year) %>% 
    
    # Group into a list of sensors by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot)) 
  
  
  
  # Format virtual sensors
  data_sensors <- db$hemiphot %>% 
    dplyr::filter(is.na(all_mapped_above) | all_mapped_above) %>% 
    dplyr::select(id_subplot, year, pacl, h_m = height) %>% 
    dplyr::left_join(
      db$subplot %>% dplyr::select(id_subplot, id_plot, x = x_m, y = y_m),
      by = "id_subplot"
    ) %>% 
    dplyr::arrange(id_subplot, year) %>% 
    dplyr::filter(year %in% years) %>% 
    dplyr::mutate(
      id_plot = paste0(prefix, "_", id_plot)
    ) %>% 
    dplyr::rename(id_sensor = id_subplot) %>% 
    
    # Group into a list of sensors by plot
    split(.$id_plot) %>% 
    purrr::map(~.x %>% dplyr:::select(-id_plot))
  
  
  # Get plot limits
  data_polygons <- db$plot %>%
    st_as_sf(wkt = "wkt") %>%
    st_geometry() %>%
    lapply(function(poly) {
      
      coords <- st_coordinates(poly)
      
      data.frame(
        x = coords[, "X"],
        y = coords[, "Y"]
      )
    })
  names(data_polygons) <- paste0(prefix, "_", db$plot$id_plot)
  
  
  # Return data IRRES
  return(
    list(
      "trees" = data_trees,
      "sensors" = data_sensors,
      "polygon" = data_polygons
    )
  )
}




# IRRES functions ----
get_irres_trees <- function(db_irres) {
  
  data_trees <- db_irres$tree_measure %>% 
    
    # Remove individuals with a single living observations in 2015
    dplyr::group_by(id_tree) %>%
    dplyr::arrange(id_tree, year) %>%
    dplyr::filter(
      !(sum(status == 0) == 1 &
          any(year == 2015 & status == 0))
    ) %>%
    dplyr::ungroup() %>%
    
    # Add multistem dbh as the mean of all stems
    dplyr::left_join(db_irres$multistem_measure %>% 
                       dplyr::filter(status == 0 | is.na(status)) %>% 
                       dplyr::group_by(id_measure) %>% 
                       dplyr::summarise(dbh_mm_multistem = mean(dbh_mm)), 
                     by = "id_measure") %>% 
    dplyr::mutate(
      dbh_mm = ifelse(multistem, dbh_mm_multistem, dbh_mm), # Multistem dbh is the mean dbh of the living stems
      dbh_mm = ifelse(year == 2015, NA, dbh_mm), # Remove 2015 DBH as they are measured with a tape
      dbh_mm = ifelse(status != 0, NA, dbh_mm), # Set non-living trees DBH as NA
      dbh_cm = dbh_mm / 10
    ) %>% 
    
    # Add plot info
    dplyr::left_join(db_irres$tree, by = "id_tree") %>% 
    dplyr::left_join(db_irres$plot, by = "id_plot") %>%
    dplyr::left_join(db_irres$dico_species, by = c("species" = "id")) %>% 
    
    dplyr::mutate(
      plot_name = paste(id_plot, plot_name, sep = "-")
    ) %>% 
    
    dplyr::select(id_plot, plot_name, x = x_m, y = y_m,
                  id_tree, speciesId = species, species_name, 
                  year, status, dbh_cm)
  
  return(data_trees)
  
} 


check_timeseries_continuity <- function(data_trees) {
  
  # 1. Define the expected previous year of the inventory time sequence
  expected_previous_years <- tidyr::crossing(
    id_plot = 1:9,
    data.frame(
      year = seq(2016, 2026, by = 2),
      expected_year_previous = c(NA, seq(2016, 2024, by = 2))
    )
  ) %>% 
    # Plot 2 has been measured in 2021 instead of 2020
    dplyr::mutate(
      year = ifelse(id_plot == 2 & year == 2020, 2021, year),
      expected_year_previous = ifelse(id_plot == 2 & year == 2022, 2021, expected_year_previous)
    )
  
  
  # 2. Check continuity of each individual tree
  continuity_check <- data_trees %>%
    dplyr::arrange(id_tree, year) %>%
    dplyr::group_by(id_tree) %>%
    dplyr::mutate(
      
      is_first_observation = row_number() == 1,
      
      year_previous = dplyr::lag(year),
      status_previous = dplyr::lag(status),
      
      # A tree is considered alive after an inventory only when
      # its previous status was 0.
      alive_previous = status_previous == 0
    ) %>% 
    dplyr::ungroup() %>% 
    
    dplyr::left_join(expected_previous_years, by = c("id_plot", "year")) %>% 
    
    # A tree has a correct time series if, for the non first observations:
    #   - The previous date is the expected one
    #   - AND the tree was alive in previous observation was alive
    dplyr::mutate(
      wrong_observation = ( is_first_observation & status != 0 ) |
        (!is_first_observation & !(year_previous == expected_year_previous & alive_previous) )
    ) %>% 
    
    # Get trees time series with one wrong observation
    dplyr::group_by(id_tree) %>% 
    dplyr::filter(sum(wrong_observation) > 0) %>% 
    dplyr::ungroup()
  
  # 3. Print for discontinuities
  # There is 5 trees with continuity issues: 
  if (any(continuity_check$wrong_observation, na.rm = TRUE)) {
    
    problematic_ids <- unique(
      continuity_check$id_tree[
        continuity_check$wrong_observation
      ]
    )
    
    message(
      "Discontinuities detected for id_tree: ",
      paste(problematic_ids, collapse = ", ")
    )
    
  } else {
    
    message("No discontinuities detected in the tree time series.")
    
  }
  
  return(data_trees)
}


clean_discontinuities <- function(data_trees) {
  
  data_trees <- data_trees %>% 
    
    #   - 373: noted as harvested in 2024, but not living in 2026, with consistent DBH measurements
    #       ==> Set the 2024 observation of harvesting as living
    dplyr::mutate(status = ifelse(id_tree == 373 & year == 2024, 0, status)) %>% 
    
    #   - 2112: the only observation is harvested tree at year 2018:
    #       ==> remove this tree
    dplyr::filter(id_tree != 2112) %>% 
    
    #   - 2670: the only observation is dead at year 2016:
    #       ==> remove this tree
    dplyr::filter(id_tree != 2670) %>% 
    
    #   - 2152: noted as dead in 2020, but living in 2024 and 2026, with consistent DBH measurements
    #       ==> Set the 2020 observation of death as living
    dplyr::mutate(status = ifelse(id_tree == 2152 & year == 2020, 0, status)) %>% 
    
    #   - 2341: noted as dead both in 2024 and in 2026, with no DBH increments, thus the tree was really dead in 2024 but I don't know why, has been measured in 2026
    #       ==> remove the 2026 record 
    dplyr::filter(!(id_tree == 2341 & year == 2026)) %>% 
    
    #   - 2841: noted as harvested in 2022, without DBH measuremenbt, thus no measurement in 2024, but a measurement in 2026 with a consistent DBH increment with 2020
    #       ==> set the tree as non harvested in 2022 and add observation for 2024: missing DBH will be estimated by individual extrapolation
    dplyr::mutate(status = ifelse(id_tree == 2841 & year == 2022, 0, status)) %>% 
    
    dplyr::add_row(id_plot = 6, plot_name = "6-Do. jeune, St-Hubert", id_tree = 2841, 
                   species_name = "Douglas", status = 0, year = 2024, dbh_cm = NA) %>% 
    
    #   - 3193: noted as dead both in 2024 and in 2026, with negative DBH increments between 2024 and 2026, thus the tree was really dead in 2024 but I don't know why, has been measured in 2026
    #       ==> remove the 2026 record
    dplyr::filter(!(id_tree == 3193 & year == 2026))
  
  
  return(data_trees)
}


estimate_dbh_2015 <- function(data_trees) {
  
  # 1. Calibration data
  data_calib <- data_trees %>%
    dplyr::filter(
      !is.na(dbh_cm),
      status == 0,
      year > 2015 # Remove data from 2015 as they are measured with a tape
    ) %>%
    dplyr::arrange(id_tree, year) %>%
    dplyr::group_by(id_tree) %>%
    dplyr::mutate(
      dbh_cm_next = dplyr::lead(dbh_cm),
      year_next = dplyr::lead(year),
      
      delta_dbh_cm = dbh_cm_next - dbh_cm,
      delta_year = year_next - year,
      
      dbh_increment_cm_year =
        delta_dbh_cm / delta_year
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(dbh_cm_next),
      !is.na(dbh_increment_cm_year),
      delta_year > 0
    )
  
  # 2. Fit DBH-increment model
  mod <- lme4::lmer(
    dbh_increment_cm_year ~
      dbh_cm_next * species_name +
      (1 | id_tree),
    data = data_calib
  )
  
  # 3. Recreate the 2015 stand
  # Use the same stand structure as in 2016 
  stand_2015 <- data_trees %>%
    dplyr::filter(year == 2015, status == 0) %>%
    dplyr::select(-dbh_cm) %>% 
    
    # Get the first living observation with non-na DBH for each tree
    dplyr::left_join(
      
      data_trees %>% 
        dplyr::filter(!is.na(dbh_cm) & status == 0) %>% 
        dplyr::group_by(id_tree) %>% 
        dplyr::arrange(id_tree, year) %>% 
        dplyr::filter(row_number() == 1) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(id_tree, year_next = year, dbh_cm_next = dbh_cm),
      by = "id_tree"
    ) %>% 
    
    
    # 4. Estimate DBH in 2015 by estimating individual yearly increment
    dplyr::mutate(
      dbh_increment_cm_yr = predict(
        mod,
        newdata = data.frame(
          dbh_cm_next = dbh_cm_next,
          id_tree = id_tree,
          species_name = species_name
        ),
        re.form = NA,
        allow.new.levels = TRUE
      ),
      
      # Avoid negative predicted increments
      dbh_increment_cm_yr = pmax(
        dbh_increment_cm_yr,
        0
      ),
      
      # Back-calculate one year
      dbh_cm_2015 = dbh_cm_next -
        dbh_increment_cm_yr / (year_next - 2015),
      
      dbh_cm_2015 = round(
        dbh_cm_2015,
        1
      )
    ) %>% 
    
    dplyr::select(id_tree, year, dbh_cm_2015)
  
  
  # 5. Add the reconstructed 2015 stands
  data_trees <- data_trees %>% 
    dplyr::left_join(stand_2015, by = c("id_tree", "year")) %>% 
    dplyr::mutate(dbh_cm = coalesce(dbh_cm, dbh_cm_2015)) %>%
    dplyr::select(-dbh_cm_2015)
  
  return(data_trees)
}


keep_living_trees <- function(data_trees) {
  
  data_trees <- data_trees %>% 
    dplyr::filter(status == 0)
  
  return(data_trees)
}


get_missing_dbh_estimate <- function(data_trees, idtree, years) {
  
  # Individual information
  dat_individual <- data_trees %>%
    dplyr::filter(
      id_tree == idtree,
      !is.na(dbh_cm),
      status == 0,
      year != 2015 # Because measured with a tape
    ) %>%
    dplyr::arrange(year)
  
  
  # Case 1: at least two individual DBH observations
  # Regression with year and all available dbh measurements of this tree
  if (nrow(dat_individual) >= 2) {
    
    mod <- lm(dbh_cm ~ year, data = dat_individual)
    
    pred <- predict(mod, newdata = data.frame(year = years))
    
    
    # Case 2: one individual DBH observation
    # Estimate DBH with mean speciesXplot dbh increment from the available DBH observation  
  } else if (nrow(dat_individual) == 1) {
    
    # Tree information
    id_plot_individual <- dat_individual$plot_name[1]
    yr_individual <- dat_individual$year[1]
    sp_individual <- dat_individual$species_name[1]
    dbh_individual <- dat_individual$dbh_cm[1]
    
    # Plot data.frame
    dat_plot <- data_trees %>%
      dplyr::filter(
        plot_name == id_plot_individual,
        !is.na(dbh_cm),
        year > 2015
      ) %>%
      dplyr::arrange(id_tree, year) %>%
      dplyr::group_by(id_tree) %>%
      dplyr::mutate(
        delta_dbh_cm = dbh_cm - dplyr::lag(dbh_cm),
        delta_year = year - dplyr::lag(year),
        
        # Annual increment in mm/year
        dbh_increment_cm_year =
          delta_dbh_cm / delta_year
      ) %>%
      dplyr::ungroup()
    
    # Mean species increment in this plot
    mean_dbh_increment_cm_year <- dat_plot %>%
      dplyr::filter(species_name == sp_individual) %>% 
      dplyr::summarise(
        mean_increment = mean(
          dbh_increment_cm_year,
          na.rm = TRUE
        )
      ) %>%
      dplyr::pull(mean_increment)
    
    # If not enough observations, keep non-species-specific value
    if (is.na(mean_dbh_increment_cm_year)) {
      mean_dbh_increment_cm_year <- dat_plot %>%
        dplyr::summarise(
          mean_increment = mean(
            dbh_increment_cm_year,
            na.rm = TRUE
          )
        ) %>%
        dplyr::pull(mean_increment)
    }
    
    # Estimate missing dbh based on the only value
    pred <- dbh_individual +
      mean_dbh_increment_cm_year *
      (years - yr_individual)
    
  } else {
    
    stop("No available individuals with id_tree = ", idtree)
  }
  
  # Return DBH in cm
  round(pred, 1)
}


estimate_missing_dbh <- function(data_trees) {
  
  # Estimate NA DBH from the dbh of other years of the same tree
  data_missingDBH <- data_trees %>% 
    dplyr::filter(is.na(dbh_cm), status == 0) %>% 
    dplyr::select(id_tree, year) %>% 
    dplyr::group_by(id_tree) %>% 
    dplyr::group_modify(~ { 
      .x %>% 
        dplyr::mutate( dbh_cm_estimate = get_missing_dbh_estimate(data_trees, .y$id_tree, year) ) 
    }) %>%
    dplyr::ungroup() 
  
  
  data_trees <- data_trees %>% 
    dplyr::left_join(data_missingDBH, by = c("id_tree", "year")) %>% 
    dplyr::mutate(
      dbh_cm = coalesce(dbh_cm, dbh_cm_estimate),
      dbh_cm = as.numeric(dbh_cm)
    ) %>% 
    dplyr::select(-dbh_cm_estimate)
  
  return(data_trees)
}


keep_adult_trees <- function(data_trees) {
  
  data_trees <- data_trees %>% 
    dplyr::filter(dbh_cm >= 10)
  
  return(data_trees)
}


keep_years <- function(data_trees, years) {
  
  data_trees <- data_trees %>% 
    dplyr::filter(year %in% years)
  
  return(data_trees)
  
}


add_crown_measurements <- function(data_trees, db_irres) {
  
  data_radius <- dplyr::bind_rows(
    get_data_radius_2026(db_irres), 
    get_data_radius_2015_2016(db_irres) # Here, crown measurements for plot 8 and 9 have mostly been done in 2015
  )
  
  data_height <- db_irres$height %>% 
    dplyr::select(id_tree, year, 
                  h_m = height_m, 
                  hbase_m = crownbase_m)
  
  
  data_trees <- data_trees %>% 
    
    dplyr::left_join(data_height, by = c("id_tree", "year")) %>% 
    dplyr::left_join(data_radius, by = c("id_tree", "year")) %>% 
    
    dplyr::mutate(
      
      rn_m = ifelse(!is.na(radius_to_treecenter_m_measured_n),
                    radius_to_treecenter_m_measured_n - dbh_cm / 100 / 2,
                    radius_m_measured_n),
      rn_m = pmax(0, rn_m),
      
      
      re_m = ifelse(!is.na(radius_to_treecenter_m_measured_e),
                    radius_to_treecenter_m_measured_e - dbh_cm / 100 / 2,
                    radius_m_measured_e),
      re_m = pmax(0, re_m),
      
      
      rs_m = ifelse(!is.na(radius_to_treecenter_m_measured_s),
                    radius_to_treecenter_m_measured_s - dbh_cm / 100 / 2,
                    radius_m_measured_s),
      rs_m = pmax(0, rs_m),
      
      
      rw_m = ifelse(!is.na(radius_to_treecenter_m_measured_w),
                    radius_to_treecenter_m_measured_w - dbh_cm / 100 / 2,
                    radius_m_measured_w),
      rw_m = pmax(0, rw_m)
    ) %>% 
    
    dplyr::select(-contains(c("radius_to_treecenter_m_measured_", "radius_m_measured")))
  
  
  return(data_trees)
}


get_data_radius_2026 <- function(db_irres) {
  
  # For 2026, we took at exactly the given azimuth
  # Each masured individual have at least four azimuth (0 = north, 90 = east, 180 = south, 270 = west)
  data_radius_2026 <- db_irres$crown %>% 
    dplyr::filter(year == 2026) %>% 
    dplyr::mutate(azimuth_str = case_match(azimuth_deg,
                                           0 ~ "n", 90 ~ "e", 180 ~ "s", 270 ~ "w"))
  
  # Id_tree 3406 need north to be estimated as there was a samll gap
  # Thus, north west and north east has been measured
  data_radius_2026 <- data_radius_2026 %>% 
    dplyr::add_row(id_tree = 3406, year = 2026, azimuth_deg = 0, azimuth_str = "n", 
                   radius_m = mean(
                     data_radius_2026$radius_m[data_radius_2026$id_tree == 3406 & data_radius_2026$azimuth_deg == 45],
                     data_radius_2026$radius_m[data_radius_2026$id_tree == 3406 & data_radius_2026$azimuth_deg == 315]
                   )) %>% 
    
    dplyr::filter(!is.na(azimuth_str)) %>% 
    
    tidyr::pivot_wider(id_cols = c(id_tree, year), names_from = azimuth_str, 
                       names_prefix = "radius_m_measured_", values_from = radius_m)
  
  return(data_radius_2026)
}


get_data_radius_2015_2016 <- function(db_irres) {
  
  data_radius_2015_2016 <- db_irres$crown %>%
    dplyr::filter(year %in% c(2015, 2016)) %>%
    dplyr::select(
      id, id_tree, year,
      xc_m = x_m, yc_m = y_m
    ) %>%
    dplyr::left_join(
      db_irres$tree %>%
        dplyr::select(id_tree, x_m, y_m),
      by = "id_tree"
    ) %>%
    dplyr::mutate(
      # Compute crown radius
      radius_to_treecenter_m = sqrt(
        (x_m - xc_m)^2 + (y_m - yc_m)^2
      ),
      
      # Compute azimuth
      azimuth_deg = (
        180 / pi * atan2(xc_m - x_m, yc_m - y_m)
      ) %% 360
    ) %>%
    
    # Assign the closest radius independently to each cardinal
    # (the same measurement can be assigned more than once)
    dplyr::select(
      id, id_tree, year,
      radius_to_treecenter_m,
      azimuth = azimuth_deg
    ) %>%
    assign_closest_cardinal() %>%

    # Format output
    dplyr::select(
      id_tree, year,
      radius_to_treecenter_m,
      azimuth_str = cardinal
    ) %>%
    tidyr::pivot_wider(
      id_cols = c(id_tree, year),
      names_from = azimuth_str,
      names_prefix = "radius_to_treecenter_m_measured_",
      values_from = radius_to_treecenter_m
    )
  
  return(data_radius_2015_2016)
}



# Global functions ----
assign_closest_cardinal <- function(df_radius) {
  
  cardinals <- c(
    n = 0,
    e = 90,
    s = 180,
    w = 270
  )
  
  df_radius %>%
    dplyr::mutate(
      azimuth = azimuth %% 360
    ) %>%
    dplyr::group_by(id_tree) %>%
    dplyr::group_modify(~ {
      
      dat <- .x
      
      # Find the closest radius independently for each cardinal
      purrr::map_dfr(names(cardinals), function(card) {
        
        # The expression always returns a difference in the range [-180, 180]
        # Without circular correction, the difference between 350 and 0 would be calculated as 350, not 10
        dist <- abs(
          ((dat$azimuth - cardinals[card] + 180) %% 360) - 180
        )
        
        best <- which.min(dist)
        
        dat[best, , drop = FALSE] %>%
          dplyr::mutate(
            cardinal = card,
            card_dist = dist[best]
          )
      })
      
    }) %>%
    dplyr::ungroup()
}

