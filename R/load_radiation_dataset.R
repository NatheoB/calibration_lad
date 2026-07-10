load_radiation_dataset <- function(fp) {
  
  df <- vroom(fp)
  
  data_rad <- split(df, df$site) %>% 
    purrr::map(~.x %>% dplyr::select(-site))
  
  data_rad
}