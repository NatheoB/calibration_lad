get_radiation_dataset <- function(data_plots) {
  
  sapply(data_plots$name, function(p) {
    
    # Get plot info
    plot_info <- data_plots %>% dplyr::filter(name == p)
    
    # Get monthly radiation from PVGIS based on latitude and longitude of the plot
    # SamsaRaLight::get_monthly_rad(
    #   longitude = plot_info$longitude,
    #   latitude = plot_info$latitude,
    # )
    SamsaRaLight::data_rad_prenovel
    
  }, USE.NAMES = TRUE, simplify = FALSE)
  
}