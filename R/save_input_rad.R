save_input_rad <- function(data_rad, fp) {
  
  df <- data_rad %>% 
    dplyr::bind_rows(.id = "site")
  
  write.table(df, fp,
              dec = ".", sep = ";", row.names = F)
  
}