clean_database <- function(init_db) {
  
  # Set plot extent without tops to NULL
  init_db$plot_extents <- init_db$plot_extents %>% 
    purrr::map_if(~nrow(.x) == 0, ~return(NULL))
  
  # Correct Baileux plot_extents
  # Set to the HULL convex whereas we want the rect
  # The bounding box in the inventories are FALSE
  init_db$plot_extents[["BaileuxBeech"]] <- data.frame(
    X = c(59.99, 66.1, -45.12, -51.23),
    Y = c(64.18, -111.27, -115.15, 60.31)
  )
  
  init_db$plot_extents[["BaileuxMixed"]] <- data.frame(
    X = c(-75.04, 84.85, 110.5, -49.39),
    Y = c(44.77, 83.64, -21.84, -60.72)
  )
  
  init_db$plot_extents[["BaileuxOak"]] <- data.frame(
    X = c(-43.82, 68.91, 55.25, -57.49),
    Y = c(67.5, 53.47, -56.32, -42.3)
  )
  
  
  # Ventoux C7, pb with species codes
  # 210 is "Autres_Feuillus" for other plots
  # But here, there are 3 occurences of 210, that are attributed to "Autre érable", "Houx" and "Genevrier"
  # Thus, set 210 species code to "Autres_Feuillus"
  init_db$species$VentouxC7 <- init_db$species$VentouxC7[which(init_db$species$VentouxC7$SpCode_SamsaraLL != 210),]
  init_db$species$VentouxC7[nrow(init_db$species$VentouxC7) + 1,] <-
    c(210, "Autre feuillu", "Autre_feuillu", "#00FF00", "VentouxC7")
  init_db$species$VentouxC7$SpCode_SamsaraLL <- as.numeric(init_db$species$VentouxC7$SpCode_SamsaraLL)
  
  
  # Set Quercus petraea as Quercus sp. in Baileux Oak
  init_db$species$BaileuxOak[which(init_db$species$BaileuxOak$SpCode_SamsaraLL == 1), 
                             c("Essence", "Essence_Latin")] <- c("Chêne indigène", "Quercus_sp")
  
  
  # Get species to calibrate and et other gymno and angio
  
  # Essence_Latin          Freq_trees
  # <chr>                 <int>
  # 1        Castanea_sativa    1
  # 2           Crataegus_sp    1
  # 3               Malus_sp    1
  # 4      Sorbus_torminalis    1
  # 5          Larix_decidua    2
  # 6             Abies_alba    3
  # 7         Autre_resineux    3
  # 8        Populus_tremula    3
  # 9       Sorbus_aucuparia    3
  # 10      Carpinus_betulus    4
  # 11             Sorbus_sp    5
  # 12 Pseudotsuga_menziesii    7
  # 13  Populus_x_canadensis    8
  # 14      Pinus_sylvestris   10
  # 15         Autre_feuillu   12
  # 16             Betula_sp   15
  # 17           Picea_abies   17
  # 18            Quercus_sp   29
  # 19       Fagus_sylvatica   36
  
  # Essence_Latin          Freq_trees
  # <chr>                 <int>
  # 1 Fagus_sylvatica        1772
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
  
  sp_gymno <- c("Picea_abies", "Pseudotsuga_menziesii", 
                "Autre_resineux", "Pinus_sylvestris",         
                "Larix_decidua", "Abies_alba")
  
  sp_angio <- c("Quercus_sp", "Fagus_sylvatica", "Autre_feuillu", 
                "Carpinus_betulus", "Betula_sp", "Populus_x_canadensis",
                "Populus_tremula", "Sorbus_sp", "Sorbus_torminalis", "Sorbus_aucuparia", 
                "Crataegus_sp", "Castanea_sativa", "Malus_sp")
    
  
  init_db$species <- init_db$species %>% 
    purrr::map(~.x %>% 
                 dplyr::mutate(
                   
                   sp_calib = case_when(
                     Essence_Latin %in% sp_to_calib ~ Essence_Latin,
                     Essence_Latin %in% sp_gymno ~ "other_gymno",
                     Essence_Latin %in% sp_angio ~ "other_angio"
                   ),
                   
                   sp_phylogeny = case_when(
                     Essence_Latin %in% sp_gymno ~ "gymnosperm",
                     Essence_Latin %in% sp_angio ~ "angiosperm"
                   )
                 )
    )
  
  
  # Remove the -1 PACL values
  init_db$sensors <- init_db$sensors %>% 
    purrr::map(~.x %>% dplyr::filter(PACLtotal <= 1, 
                                     PACLtotal >= 0))
  
  
  # Change some wrong trees info
  
  # ---- Creating SamsaRaLight stand for Lorris38..
  # `hbase_m` must be strictly lower than `h_m` for tree(s) with id_tree: 87
  lorris38_pisy <- subset(init_db[["trees"]][["Lorris38"]], SpCode == 45 & Id != 87)
  mod_dbh_h_pisy <- lm(H ~ Dbh, lorris38_pisy)
  # summary(mod_dbh_h_pisy)
  
  dbh_lorris38_id87 <- subset(init_db[["trees"]][["Lorris38"]], Id == 87, select = Dbh)
  h_lorris38_id87 <- predict(mod_dbh_h_pisy, dbh_lorris38_id87)
  
  init_db[["trees"]][["Lorris38"]]$H[which(init_db[["trees"]][["Lorris38"]]$Id == 87)] <- h_lorris38_id87
  
  
  # Return final database
  return(init_db)
}