load_database_from_rdata <- function(rdata_fp) {
  
  # Load data into a controlled R environment
  e <- new.env()
  load(rdata_fp, envir = e)
  
  # Create plots dataset
  data_plots <- e$SensorsDf %>% 
    dplyr::select(name = plot, origin = origine) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(
      longitude = 5.2,
      latitude = 50.04
    ) %>% 
    as.data.frame()
  

  # Get data and save lists into a list of datasets
  # For each list as a list of sites, rename each element by the site name
  db_list <- list(
    "plots" = data_plots,
    "plot_extents" = setNames(e$listPlotExtent, data_plots$name),
    "sensors" = setNames(e$listSensors, data_plots$name),
    "trees" = setNames(e$listTree, data_plots$name),
    "species" = setNames(e$listSp, data_plots$name)
  )
  
  
  # Set plot extent without tops to NULL
  db_list$plot_extents <- db_list$plot_extents %>% 
    purrr::map_if(~nrow(.x) == 0, ~return(NULL))
  
  # Ventoux C7, pb with species codes
  # 210 is "Autres_Feuillus" for other plots
  # But here, there are 3 occurences of 210, that are attributed to "Autre érable", "Houx" and "Genevrier"
  # Thus, set 210 species code to "Autres_Feuillus"
  db_list$species$VentouxC7 <- db_list$species$VentouxC7[which(db_list$species$VentouxC7$SpCode_SamsaraLL != 210),]
  db_list$species$VentouxC7[nrow(db_list$species$VentouxC7) + 1,] <-
    c(210, "Autre feuillu", "Autre_feuillu", "#00FF00", "VentouxC7")
  db_list$species$VentouxC7$SpCode_SamsaraLL <- as.numeric(db_list$species$VentouxC7$SpCode_SamsaraLL)
  
  
  # Return final database
  return(db_list)
}