get_radiation_dataset <- function(data_plots) {
  
  sapply(data_plots$name, function(p) {
    
    message("--- Fetching radiations for ", p, "...")
    
    # Get plot info
    plot_info <- data_plots %>% dplyr::filter(name == p)
    
    # Get monthly radiation from PVGIS based on latitude and longitude of the plot
    rad <- SamsaRaLight::get_monthly_radiations(
      longitude = plot_info$longitude,
      latitude = plot_info$latitude,
    )

    # Test validity of the radiation dataset
    SamsaRaLight::check_monthly_radiations(rad)
    
    # Return the radiation dataset
    rad
    
  }, USE.NAMES = TRUE, simplify = FALSE)
  
}