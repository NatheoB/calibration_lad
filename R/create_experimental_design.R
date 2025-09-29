create_experimental_design <- function(site_names, n_replicate, lad_values) {
  
  expand.grid(
    site = site_names,
    replicate = 1:n_replicate
  ) %>% 
    tidyr::unite(
      col = id,
      site, replicate,
      sep = "_",
      remove = FALSE
    ) %>% 
    dplyr::relocate(id)
  
}
