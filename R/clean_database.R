clean_database <- function(init_db_raw) {
  
  # Set plot extent without tops to NULL
  init_db_raw$plot_extents <- init_db_raw$plot_extents %>% 
    purrr::map_if(~nrow(.x) == 0, ~return(NULL))
  
  # Correct Baileux plot_extents
  # Set to the HULL convex whereas we want the rect
  # The bounding box in the inventories are FALSE
  init_db_raw$plot_extents[["BaileuxBeech"]] <- data.frame(
    X = c(59.99, 66.1, -45.12, -51.23),
    Y = c(64.18, -111.27, -115.15, 60.31)
  )
  
  init_db_raw$plot_extents[["BaileuxMixed"]] <- data.frame(
    X = c(-75.04, 84.85, 110.5, -49.39),
    Y = c(44.77, 83.64, -21.84, -60.72)
  )
  
  init_db_raw$plot_extents[["BaileuxOak"]] <- data.frame(
    X = c(-43.82, 68.91, 55.25, -57.49),
    Y = c(67.5, 53.47, -56.32, -42.3)
  )
  
  
  # Ventoux C7, pb with species codes
  # 210 is "Autres_Feuillus" for other plots
  # But here, there are 3 occurences of 210, that are attributed to "Autre érable", "Houx" and "Genevrier"
  # Thus, set 210 species code to "Autres_Feuillus"
  init_db_raw$species$VentouxC7 <- init_db_raw$species$VentouxC7[which(init_db_raw$species$VentouxC7$SpCode_SamsaraLL != 210),]
  init_db_raw$species$VentouxC7[nrow(init_db_raw$species$VentouxC7) + 1,] <-
    c(210, "Autre feuillu", "Autre_feuillu", "#00FF00", "VentouxC7")
  init_db_raw$species$VentouxC7$SpCode_SamsaraLL <- as.numeric(init_db_raw$species$VentouxC7$SpCode_SamsaraLL)
  
  
  # Set Quercus petraea as Quercus sp. in Baileux Oak
  init_db_raw$species$BaileuxOak[which(init_db_raw$species$BaileuxOak$SpCode_SamsaraLL == 1), 
                             c("Essence", "Essence_Latin")] <- c("Chêne indigène", "Quercus_sp")
  
  
  # Get species in the inventory
  message("Species in unique stands...")
  
  print(
    init_db_raw$species %>%
      purrr::map(~unique(.x$Essence_Latin)) %>%
      unlist() %>%
      unname() %>%
      table(useNA = "always") %>%
      sort()
  )
  
  # Castanea_sativa          Crataegus_sp              Malus_sp       Quercus_petraea     Sorbus_torminalis 
  # 1                     1                     1                     1                     1 
  # Larix_decidua            Abies_alba        Autre_resineux       Populus_tremula      Sorbus_aucuparia 
  # 2                     3                     3                     3                     3 
  # Carpinus_betulus             Sorbus_sp Pseudotsuga_menziesii  Populus_x_canadensis      Pinus_sylvestris 
  # 4                     5                     7                     8                    10 
  # Autre_feuillu             Betula_sp           Picea_abies            Quercus_sp       Fagus_sylvatica 
  # 12                    15                    17                    28                    36 
  
  message("Number of trees of each species...")
  
  print(
    init_db_raw$trees %>% 
      purrr::map2(names(init_db_raw$trees), 
                  ~.x %>% 
                    dplyr::left_join(init_db_raw$species[[.y]], 
                                     by = c("SpCode" = "SpCode_SamsaraLL")) %>% 
                    dplyr::pull(Essence_Latin)) %>% 
      unlist() %>%
      unname() %>%
      table(useNA = "always") %>%
      sort()
  )
  
  # Castanea_sativa          Crataegus_sp     Sorbus_torminalis              Malus_sp        Autre_resineux 
  # 1                     1                     1                     5                    13 
  # Populus_tremula             Sorbus_sp         Larix_decidua         Autre_feuillu            Abies_alba 
  # 23                    39                    73                    78                   102 
  # Betula_sp  Populus_x_canadensis      Sorbus_aucuparia       Quercus_petraea      Carpinus_betulus 
  # 112                   112                   132                   189                   288 
  # Pinus_sylvestris Pseudotsuga_menziesii            Quercus_sp           Picea_abies       Fagus_sylvatica 
  # 866                  1094                  1497                  1771                  1772 
  
  
  # Specify species to fot and wheter the species is angio or gymnosperm
  sp_to_calib <- list(
    "Fagus_sylvatica" = "Fagus_sylvatica", 
    "Picea_abies" = "Picea_abies", 
    "Quercus_sp" = c("Quercus_sp", "Quercus_petraea"), 
    "Pseudotsuga_menziesii" = "Pseudotsuga_menziesii",
    "Pinus_sylvestris" = "Pinus_sylvestris", 
    "Carpinus_betulus" = "Carpinus_betulus", 
    "Sorbus_sp" = c("Sorbus_aucuparia", "Sorbus_sp", "Sorbus_torminalis"),
    "Populus_sp" = c("Populus_x_canadensis", "Populus_tremula"),
    "Betula_sp" = "Betula_sp",
    "Abies_alba" = "Abies_alba", 
    "Larix_decidua" = "Larix_decidua",
    "other_angio" = c("Autre_feuillu", "Malus_sp", "Crataegus_sp", "Castanea_sativa"),
    "other_gymno" = "Autre_resineux"
  )
  sp_to_calib_df <- tibble(
    sp_calib = names(sp_to_calib),
    species_in_inv = sp_to_calib
  ) %>%
    tidyr::unnest(species_in_inv)
  
  sp_gymno <- c(
    "Picea_abies", 
    "Pseudotsuga_menziesii", 
    "Pinus_sylvestris",         
    "Abies_alba",
    "Larix_decidua", 
    "Autre_resineux"
  )
  
  sp_angio <- c(
    "Fagus_sylvatica", 
    "Quercus_sp",
    "Quercus_petraea",
    "Carpinus_betulus", 
    "Sorbus_sp", 
    "Sorbus_aucuparia", 
    "Sorbus_torminalis", 
    "Populus_x_canadensis",
    "Populus_tremula", 
    "Betula_sp", 
    "Autre_feuillu",
    "Malus_sp",
    "Crataegus_sp", 
    "Castanea_sativa"
  )
  
  
  init_db_raw$species <- init_db_raw$species %>% 
    purrr::map(~.x %>% 
                 dplyr::left_join(
                   sp_to_calib_df, by = c("Essence_Latin" = "species_in_inv")
                 ) %>% 
                 dplyr::mutate(
                   functional_group = case_when(
                     Essence_Latin %in% sp_gymno ~ "gymnosperm",
                     Essence_Latin %in% sp_angio ~ "angiosperm"
                   )
                 )
    )
  
  
  # Remove the -1 PACL values
  init_db_raw$sensors <- init_db_raw$sensors %>% 
    purrr::map(~.x %>% dplyr::filter(PACLtotal <= 1, 
                                     PACLtotal >= 0))
  
  
  # Change some wrong trees info
  
  # ---- Creating SamsaRaLight stand for Lorris38..
  # `hbase_m` must be strictly lower than `h_m` for tree(s) with id_tree: 87
  lorris38_pisy <- subset(init_db_raw[["trees"]][["Lorris38"]], SpCode == 45 & Id != 87)
  mod_dbh_h_pisy <- lm(H ~ Dbh, lorris38_pisy)
  # summary(mod_dbh_h_pisy)
  
  dbh_lorris38_id87 <- subset(init_db_raw[["trees"]][["Lorris38"]], Id == 87, select = Dbh)
  h_lorris38_id87 <- predict(mod_dbh_h_pisy, dbh_lorris38_id87)
  
  init_db_raw[["trees"]][["Lorris38"]]$H[which(init_db_raw[["trees"]][["Lorris38"]]$Id == 87)] <- h_lorris38_id87
  
  
  # Return final database
  return(init_db_raw)
}