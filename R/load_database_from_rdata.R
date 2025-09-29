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
  
  
  # Set Quercus petraea as Quercus sp. in Baileux Oak
  db_list$species$BaileuxOak[which(db_list$species$BaileuxOak$SpCode_SamsaraLL == 1), c("Essence", "Essence_Latin")] <- c("Chêne indigène", "Quercus_sp")
  
  # Get species to calibrate
  
  # # A tibble: 19 × 2
  # Essence_Latin          Freq
  # <chr>                 <int>
  #   1 Fagus_sylvatica        1772
  # 2 Picea_abies            1771
  # 3 Quercus_sp             1686
  # 4 Pseudotsuga_menziesii  1094
  # 5 Pinus_sylvestris        866
  # 6 Carpinus_betulus        288
  # 7 Sorbus_aucuparia        132
  # 8 Betula_sp               112
  # 9 Populus_x_canadensis    112
  # 10 Abies_alba              102
  # 11 Autre_feuillu            78
  # 12 Larix_decidua            73
  # 13 Sorbus_sp                39
  # 14 Populus_tremula          23
  # 15 Autre_resineux           13
  # 16 Malus_sp                  5
  # 17 Castanea_sativa           1
  # 18 Crataegus_sp              1
  # 19 Sorbus_torminalis         1
  
  # Abies alba and Larix decidua because all the individuals are concentrate within few plots
  # So we will try to fit them
  sp_to_calib <- c("Fagus_sylvatica", "Picea_abies", "Quercus_sp", "Pseudotsuga_menziesii",
                   "Pinus_sylvestris", "Carpinus_betulus", "Abies_alba", "Larix_decidua")
  
  db_list$species <- db_list$species %>% 
    purrr::map(~.x %>% 
                 dplyr::mutate(sp_calib = case_when(
                   Essence_Latin %in% sp_to_calib ~ Essence_Latin,
                   TRUE ~ "other"
                 )))
  
  # Return final database
  return(db_list)
}
