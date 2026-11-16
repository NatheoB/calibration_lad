plot_stands <- function(data_stands,
                        species_colors_fp,
                        output_plots_fp) {
  
  # Set colors
  species_colors <- vroom::vroom(species_colors_fp, show_col_types = F)
  
  species_colors_inv <- setNames(
    species_colors$color_inv, 
    species_colors$sp_in_inv
  )
  
  species_colors_calib <- setNames(
    unique(species_colors$color_calib), 
    unique(species_colors$sp_calib)
  )
  
  # Create output folder
  if (!file.exists(output_plots_fp)) {
    dir.create(output_plots_fp, recursive = T) 
  }
  
  
  
  # Create the plots of each stand  
  for (site in names(data_stands)) {
    
    message("---- Plotting SamsaRaLight stand for ", site, "...")
    
    # Get the stand
    tmp_stand <- data_stands[[site]]
  
    
    # Plot stand map
    message("SamsaRaLight stand map")

    png(file.path(output_plots_fp,
                  paste0(site, "_map_sl.png")),
        width = 600, height = 600)

    print(
      plot(tmp_stand) &
        scale_fill_manual(values = species_colors_calib)
    )

    dev.off()
    
    
    # Plot stand top-down
    message("Top-down SamsaRaLight stand")
    
    png(file.path(output_plots_fp,
                  paste0(site, "_topdown_sl.png")), 
        width = 600, height = 400)
    
    print(
      plot(tmp_stand, top_down = TRUE) & 
        scale_color_manual(values = species_colors_calib)
    )
    
    dev.off()
    
    
    ## Plot inventory map
    message("Inventory map")
    
    # Change species name for plotting inventoried species
    tmp_stand$inventory <- tmp_stand$inventory %>%
      dplyr::mutate(species = sp_in_inv)
    
    
    png(file.path(output_plots_fp,
                  paste0(site, "_inv.png")), 
        width = 600, height = 600)
    
    print(
      plot_inventory(tmp_stand$inventory) + 
        scale_fill_manual(values = species_colors_inv) +
        labs(title = paste("Inventory map of", site)) +
        theme(plot.title = element_text(hjust = 0.5))
    )
    
    dev.off()
    
  }
  
  output_plots_fp
}