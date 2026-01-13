fit_lad_method1 <- function(output_sensors,
                            convergence_threshold,
                            data_stands,
                            output_punobs,
                            output_plots_fp) 
{
  
  # Erase and create folder
  unlink(output_plots_fp, recursive = T)
  dir.create(output_plots_fp, recursive = T) # create folder
  
  # Observed pacl for each site
  data_sensors <- data_stands %>% 
    purrr::map(~.x$sensors) %>% 
    dplyr::bind_rows(.id = "site") %>% 
    dplyr::left_join(output_punobs, by = c("site", "id_sensor"))
  
  # Compute residuals for each sensor and LAD
  output_residuals <- output_sensors %>% 
    dplyr::left_join(data_sensors, 
                     by = c("site", "id_sensor")) %>% 
    dplyr::mutate(
      res = PACLtotal - pacl
    )
  
  info_sensors <- output_residuals %>% 
    dplyr::group_by(site, id_sensor) %>% 
    dplyr::summarise(
      best_lad = lad[which.min(abs(res))]
    ) %>% 
    dplyr::mutate(
      converged = best_lad < convergence_threshold
    ) %>% 
    dplyr::left_join(data_sensors, 
                     by = c("site", "id_sensor")) %>% 
    dplyr::mutate(
      sensor_label = paste0(id_sensor, 
                            " - PACL = ", round(PACLtotal, 2), 
                            " - punobs = ", round(punobs, 2))
    )
  
  output_residuals_labeled <- output_residuals %>% 
    dplyr::left_join(info_sensors %>% 
                       dplyr::select(site, id_sensor, sensor_label),
                     by = c("site", "id_sensor"))
  
  # Save graph for each site
  site_names <- unique(output_sensors$site)
  
  pb <- txtProgressBar(min = 0, max = length(site_names), 
                       style = 3, width = 50, char = "=")
  i_pb <- 0
  
  for (site_name in site_names) {
    
    ## Get the site outputs
    tmp_out_site <- output_residuals_labeled %>% 
      dplyr::filter(site == site_name)
    
    tmp_info_sensors <- info_sensors %>% 
      dplyr::filter(site == site_name)
    
    ## Get optimum LAD
    optim_lad <- mean(tmp_info_sensors %>% 
                        dplyr::filter(converged) %>% 
                        dplyr::pull(best_lad))
    
    ## Set colors
    colors_converged <- c("salmon", "forestgreen")
    if (sum(tmp_info_sensors$converged) == nrow(tmp_info_sensors)) {
      colors_converged <- c("forestgreen")
    }
    
    ## Plot the graph
    plt_site_res <- ggplot() +
      facet_wrap(~sensor_label) +
      geom_line(data = tmp_out_site,
                mapping = aes(y = res, x = lad)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(data = tmp_info_sensors,
                 mapping = aes(xintercept = best_lad,
                               color = converged),
                 linewidth = 1.05) +
      geom_label(data = tmp_info_sensors,
                 mapping = aes(label = best_lad,
                               color = converged),
                 x = 3.5, y = -0.4,
                 size = 6) +
      scale_color_manual(values = colors_converged) +
      ylab("Residuals between observed and predicted total PACL") +
      xlab("Site-specific mean species LAD") +
      labs(title = paste("LAD calibration for", site_name),
           subtitle = paste0(sum(tmp_info_sensors$converged), "/", nrow(tmp_info_sensors),
                             " sensors have converged - optimum LAD = ",
                             round(optim_lad, 3))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 30),
            plot.subtitle = element_text(hjust = 0.5, size = 24))
    
    ## Save into the folder
    filepath <- file.path(output_plots_fp,
                          paste0("output_sitespecificLAD_", site_name, ".png"))
    
    png(filepath, width = 1600, height = 1200)
    print(plt_site_res)
    dev.off()
    
    ## Update pb
    i_pb <- i_pb + 1
    setTxtProgressBar(pb, i_pb)
  }
  close(pb)
  
  # Return the summarized info per sensor
  info_sensors %>% 
    dplyr::select(site, id_sensor, best_lad, converged, punobs, pacl_obs = PACLtotal)
}