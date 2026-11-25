create_sl_geometries <- function(data_stands,
                                 data_rad) {
  
  # Output list
  out_geometries <- setNames(
    vector("list", length = length(data_stands)),
    names(data_stands)
  )
  
  # For each site, precompute the SamsaRaLight geometry (rays x crowns, sensors only)
  # Used in the likelihood to estimate PACL for any set of tree LAD (see SamsaRaLight::estimate_pacl_from_geometry())
  for (site in names(out_geometries)) {
    
    message("---- Precomputing SamsaRaLight geometry for ", site, "...")
    
    out_geometries[[site]] <- SamsaRaLight::precompute_sl_geometry(
      data_stands[[site]],
      data_rad[[site]],
      sensors_only = TRUE,
      parallel_mode = TRUE,
      verbose = FALSE
    )
  }
  
  out_geometries
}
