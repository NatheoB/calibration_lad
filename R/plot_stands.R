plot_stands <- function(data_stands,
                        data_species,
                        output_plots_fp) {
  
  
  # Create the plots of each stand  
  for (site in names(data_stands)) {
    
    message("---- Plotting SamsaRaLight stand for ", site, "...")
    
    # Get the stand
    tmp_stand <- data_stands[[site]]
    
    
    # Set colors
    tmp_colors_inv <- setNames(data_species[[site]]$Couleur, 
                               data_species[[site]]$Essence_Latin)
    
    tmp_colors_stand <- setNames(data_species[[site]]$Couleur, 
                                 data_species[[site]]$sp_calib)
    tmp_colors_stand[c("other_angio", "other_gymno")] <- c("#00FF00", "#D2691E")
    
    
    # Create output folder
    if (!file.exists(output_plots_fp)) {
      dir.create(output_plots_fp, recursive = T) 
    }
    
    
    # Plot stand map
    message("SamsaRaLight stand map")
    
    png(file.path(output_plots_fp,
                  paste0(site, "_map_sl.png")), 
        width = 600, height = 600)
    
    print(
      plot(tmp_stand) & 
        scale_fill_manual(values = tmp_colors_stand)
    )
    
    dev.off()
    
    
    # Plot stand top-down
    message("Top-down SamsaRaLight stand")
    
    png(file.path(output_plots_fp,
                  paste0(site, "_topdown_sl.png")), 
        width = 600, height = 400)
    
    print(
      plot(tmp_stand, top_down = TRUE) & 
        scale_color_manual(values = tmp_colors_stand)
    )
    
    dev.off()
    
    
    ## Plot inventory map
    message("Inventory map")
    
    # Change species name for plotting inventoried species
    tmp_stand$inventory <- tmp_stand$inventory %>%
      dplyr::mutate(species = species_in_inv)
    
    
    png(file.path(output_plots_fp,
                  paste0(site, "_inv.png")), 
        width = 600, height = 600)
    
    print(
      plot_inventory(tmp_stand$inventory) + 
        scale_fill_manual(values = tmp_colors_inv) +
        labs(title = paste("Inventory map of", site)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
    
    dev.off()
    
  }
  
  output_plots_fp
}