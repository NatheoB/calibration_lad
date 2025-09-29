plot_stand <- function(data_stand, data_species,
                       output_folderpath, plot_name) {
  
  # Initialize datasets ----
  
  ## Add species information to the tree dataset
  data_trees <- data_stand$trees %>% 
    dplyr::left_join(data_species,
                     by = c("species_code" = "SpCode_SamsaraLL")) %>% 
    dplyr::mutate(Essence_Latin = factor(Essence_Latin, levels = data_species$Essence_Latin))
  
  
  ## Create cell dataframe
  cells_xcenter <- seq(data_stand$info$cell_size / 2, 
                       data_stand$info$cell_size * data_stand$info$n_cells_x - data_stand$info$cell_size / 2,
                       by = data_stand$info$cell_size)
  
  cells_ycenter <- seq(data_stand$info$cell_size / 2, 
                       data_stand$info$cell_size * data_stand$info$n_cells_y - data_stand$info$cell_size / 2,
                       by = data_stand$info$cell_size)
  
  data_cells <- expand.grid(
    x_center = cells_xcenter,
    y_center = cells_ycenter
  )
  
  
  # Create the plot ----
  stand_plot <- ggplot() +
    coord_equal() +
    
    # CELLS 
    geom_tile(data = data_cells, 
              mapping = aes(x = x_center, y = y_center), 
              fill = "white", color = "darkgray")

  
    # CORE POLYGON
  if (!is.null(data_stand$core_polygon)) {  
    stand_plot <- stand_plot +
      geom_polygon(data = data_stand$core_polygon,
                   mapping = aes(x = x, y = y),
                   fill = "yellow", color = "black", alpha = 0.5)
  }

    
    # TREES
  stand_plot <- stand_plot +
    geom_ellipse(data = data_trees, 
                 mapping = aes(x0 = x, y0 = y, 
                               a = (re_m + rw_m) / 2,
                               b = (rn_m + rs_m) / 2,
                               angle = 0,
                               fill = Essence_Latin)) +
    scale_fill_manual(values = data_species$Couleur) +
    labs(fill = "Species") +
    
    # SENSORS
    geom_rect(data = data_stand$sensors,
              mapping = aes(xmin = x - 0.5,
                            ymin = y - 0.5,
                            xmax = x + 0.5,
                            ymax = y + 0.5),
              color = "black", fill = "black") +
    
    # GRAPHIC
    scale_x_continuous(breaks = seq(0, data_stand$info$n_cells_x * data_stand$info$cell_size,
                                    by = data_stand$info$cell_size),
                       labels = round(seq(0, data_stand$info$n_cells_x * data_stand$info$cell_size,
                                          by = data_stand$info$cell_size),
                                      digits = 1)) +
    scale_y_continuous(breaks = seq(0, data_stand$info$n_cells_y * data_stand$info$cell_size,
                                    by = data_stand$info$cell_size),
                       labels = round(seq(0, data_stand$info$n_cells_y * data_stand$info$cell_size,
                                          by = data_stand$info$cell_size),
                                      digits = 1)) +
    xlab("") + ylab("") +
    
    labs(title = plot_name,
         subtitle = paste0(round(data_stand$info$core_area_ha, 2), "ha - ",
                           round(data_stand$info$core_batot_m2ha, 2), "m2/ha - ",
                           nrow(data_stand$trees), "trees")) +
    
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) 
  
  
  # Save the plot ----
  
  ## Create folder if it does not exist
  if (!file.exists(output_folderpath)) {
    dir.create(output_folderpath, recursive = T) # create folder
  }
  
  # Save into the folder
  filepath <- file.path(output_folderpath,
                        paste0(plot_name, ".png"))
  
  png(filepath, width = 600, height = 600)
  print(stand_plot)
  dev.off()
  
}
