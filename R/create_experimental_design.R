create_experimental_design <- function(inv_names, lads, species) {
  
  data.frame(
    inv_id = paste0("inv_", rep(1:length(lads), each = length(species))),
    inv_name = inv_names,
    sp = rep(species, times = length(lads)),
    lad = rep(lads, each = length(species))
  )
  
}
