get_occurences_species2calib <- function(init_db) {
  
  init_db$trees %>% 
    purrr::map2(names(init_db$trees),
                ~.x %>% 
                  dplyr::left_join(init_db$species[[.y]],
                                   by = c("SpCode"= "SpCode_SamsaraLL"))) %>% 
    purrr::map(~table(.x$sp_calib) %>% 
                 as.data.frame() %>% 
                 rename(species = Var1)) %>% 
    dplyr::bind_rows(.id = "site")
  
}