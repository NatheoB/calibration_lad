get_sensors_punobs <- function(init_db, 
                               data_calib,
                               data_rad,
                               thresholds_punobs,
                               output_plots_fp) {
  
  # Recreate folder
  unlink(output_plots_fp, recursive = T)
  dir.create(output_plots_fp, recursive = T)
  
  
  # Apply SamsaRaLight on the sites
  site_names <- init_db$plots$name
  
  # For each site, run SamsaraLight and plot sensors interceptions ----
  pb <- txtProgressBar(min = 0, max = length(site_names), 
                       style = 3, width = 50, char = "=")
  i_pb <- 0
  
  out_sl <- setNames(vector("list", length(site_names)), site_names)
  for (site in site_names) {
    
    # Run SamsaraLight
    tmp_plot <- init_db$plots %>% 
      dplyr::filter(name == site)
    
    out_sl[[site]] <- 
      SamsaRaLight::sl_run(data_calib[[site]]$trees, 
                           data_rad[[site]],
                           sensors = data_calib[[site]]$sensors, 
                           sensors_only = FALSE,
                           latitude = tmp_plot$latitude, 
                           slope = tmp_plot$slope, 
                           aspect = tmp_plot$aspect, 
                           north_to_x_cw = tmp_plot$northToX,
                           start_day = 121, 
                           end_day = 273,
                           cell_size = data_calib[[site]]$info$cell_size, 
                           n_cells_x = data_calib[[site]]$info$n_cells_x, 
                           n_cells_y = data_calib[[site]]$info$n_cells_y,
                           turbid_medium = TRUE,
                           trunk_interception = TRUE,
                           soc = TRUE,
                           height_anglemin = 15,
                           direct_startoffset = 0,
                           direct_anglestep = 5,
                           diffuse_anglestep = 15,
                           detailed_output = TRUE
      )
    
    # Get number of interceptions for each sensor and tree ----
    interceptions_sensors <- out_sl[[site]]$output$interceptions$sensors$number$total

    # # Plot interceptions of each tree around a sensor ----
    # id_sensors <- colnames(interceptions_sensors)
    # for (id in id_sensors) {
    # 
    #   tmp_output <- out_sl[[site]]
    # 
    #   ## Bind tree dataset and interceptions ----
    #   tmp_output$output$trees <- tmp_output$output$trees %>%
    #     dplyr::left_join(
    #       data.frame(
    #         id_tree = as.integer(rownames(interceptions_sensors)),
    #         n_interceptions = unname(interceptions_sensors[,id]),
    #         prop_interceptions = unname(interceptions_sensors[,id]) / sum(interceptions_sensors[,id])
    #       ),
    #       by = "id_tree"
    #     )
    # 
    #   ## Filter the given sensor ----
    #   tmp_output$output$sensors <- tmp_output$output$sensors %>%
    #     dplyr::filter(id_sensor == as.integer(id))
    # 
    # 
    #   ## Plot the output light and interceptions ----
    #   tmp_plt_sensor <- plot_sl_output(
    #     sl_output = tmp_output,
    #     trees.fill = "n_interceptions",
    #     trees.fill.palette = "viridis",
    #     cells.fill = "pacl_slope",
    #     cells.fill.palette = "light",
    #     cells.fill.limits = c(0, 1),
    #     sensors.plot = TRUE
    #   ) +
    #     labs(title = paste0(site, " - sensor n°", id),
    #          subtitle = paste0("Part of total PACL unobstructed = ",
    #                            round(tmp_output$output$sensors$punobs_horizontal * 100, 1), "%")) +
    #     theme(plot.title = element_text(hjust = 0.5),
    #           plot.subtitle = element_text(hjust = 0.5))
    # 
    #   ## Save into the folder
    #   filepath <- file.path(output_plots_fp,
    #                         paste0(site, "_sensor", id, ".png"))
    # 
    #   png(filepath, width = 800, height = 400)
    #   print(tmp_plt_sensor)
    #   dev.off()
    # }
    
    i_pb <- i_pb + 1
    setTxtProgressBar(pb, i_pb)
  }
  close(pb)
  
  
  # Bind sensors output ----
  out_sl_sensors <- dplyr::bind_rows(
    purrr::map(out_sl, ~.x$output$sensors),
    .id = "site"
  ) %>% 
    dplyr::left_join(
      init_db$sensors %>% dplyr::bind_rows() %>% 
        dplyr::select(site = plot, id_sensor = id, 
                      PACLobs = PACLtotal,
                      PACLobs_direct = PACLdirect,
                      PACLobs_diffuse = PACLdiffuse),
      by = c("site", "id_sensor")
    ) %>% 
    dplyr::select(
      site, id_sensor, 
      PACLobs, PACLobs_direct, PACLobs_diffuse,
      punobs_horizontal, punobs_horizontal_direct, punobs_horizontal_diffuse
    )
  
  
  # For each thresholds of unobstructed PACL, plot and get the sensors to keep 
  out_thresholds <- setNames(vector("list", length(thresholds_punobs)), as.character(thresholds_punobs))
  for (trsh in thresholds_punobs) {
    
    # Get the sensors to kept for each site
    data_sensors_tokeep <- out_sl_sensors %>% 
      dplyr::mutate(to_keep = punobs_horizontal < trsh) %>%  
      dplyr::group_by(site) %>% 
      dplyr::mutate(site_keep = paste0(site, " (", sum(to_keep), "/", n(), ")")) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(above_threshold = if_else(to_keep, 
                                              paste0("Below the ", trsh * 100, "% threshold on PACL unobstructed"),
                                              paste0("Above the ", trsh * 100, "% threshold on PACL unobstructed")))
    
    # Plot the filter applied
    
    ## Create the plot
    plot_tokeep <- ggplot(data_sensors_tokeep, 
                          aes(y = punobs_horizontal * 100, x = PACLobs * 100, 
                              color = site_keep,
                              shape = above_threshold)) +
      geom_point() +
      scale_shape_manual(values = c(1, 16)) +
      geom_hline(yintercept = trsh * 100, linetype = "dashed") +
      xlab("Observed total PACL (%)") +
      ylab("Part of total PACL unobstructed (%)") +
      labs(shape = "",  
           color = "Sites (number of sensors kept/total)",
           title = paste0("Keep sensors with less than ", trsh * 100, "% of unobstructed to observed PACL"),
           subtitle = paste0(sum(data_sensors_tokeep$to_keep), "/", nrow(data_sensors_tokeep), " sensors kept")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    ## Save into the folder
    filepath <- file.path(output_plots_fp,
                          paste0("_sensors_filtering_", trsh*100, ".png"))
    
    png(filepath, width = 800, height = 400)
    print(plot_tokeep)
    dev.off()
  }
  
  # Return the sensors with estimation of PACL unobstructed
  out_sl_sensors
}