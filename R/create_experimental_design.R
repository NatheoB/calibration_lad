create_experimental_design <- function(inv_names, n_replicate, lad_values) {
  
  expand.grid(
    inv_name = inv_names,
    replicate = 1:n_replicate,
    lad = lad_values
  ) %>% 
    tidyr::unite(
      col = id,
      inv_name, replicate, lad,
      sep = "_",
      remove = FALSE
    ) %>% 
    dplyr::relocate(id)
  
}
