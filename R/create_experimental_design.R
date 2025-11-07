create_experimental_design <- function() {
  
  data.frame(
    id_mod = 1,
    species_specific = c(FALSE),
    site_rd_effect = c(FALSE),
    origin_rd_effect = c(FALSE),
    dbh_effect = c(FALSE),
    compet_effect = c(FALSE)
  )
  
  # Mettre compet variable aussi
  
  # expand.grid(
  #   species_specific = c(FALSE, TRUE),
  #   dbh_effect = c(FALSE, TRUE),
  #   compet_effect = c(FALSE, TRUE)
  # ) %>%
  #   dplyr::mutate(id_mod = row_number()) %>%
  #   dplyr::relocate(id_mod)
  
}
