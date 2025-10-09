compute_competition_trees <- function(data_trees, plot_area) {
  
  # Compute the weight to apply to each tree's basal area
  # Allowing to standardize basal area per hectare
  plot_weight <- 1 / plot_area
  
  # Compute competition variables
  data_trees %>% 
    
    # Compute basal area per hectare of each tree
    # Convert DBH (centimeters) in meters
    dplyr:::mutate(G_ha = (pi * (dbh_cm/100)^2 / 4) * plot_weight) %>% 
  
    # Compute BAtot (total basal area of the stand)
    dplyr::mutate(batot_m2ha = sum(G_ha)) %>% 
  
    # Compute BALd (basal area of bigger diameters) 
    dplyr::arrange(desc(dbh_cm), .by_group = TRUE) %>% 
    dplyr::mutate(bald_m2ha = ave(G_ha, FUN = cumsum) - G_ha) %>% 
    
    # Compute BALh (basal area of bigger heights) 
    dplyr::arrange(desc(h_m), .by_group = TRUE) %>% 
    dplyr::mutate(balh_m2ha = ave(G_ha, FUN = cumsum) - G_ha)
  
}