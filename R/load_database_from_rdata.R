load_database_from_rdata <- function(init_db_fp,
                                     plot_infos_fp) {
  
  # Get plot dataset
  data_plots <- vroom(plot_infos_fp, show_col_types = FALSE)
  
  # Load data into a controlled R environment
  e <- new.env()
  load(init_db_fp, envir = e)
  
  # Get data and save lists into a list of datasets
  # For each list as a list of sites, rename each element by the site name
  db_list <- list(
    "plots" = data_plots,
    "plot_extents" = setNames(e$listPlotExtent, data_plots$name),
    "sensors" = setNames(e$listSensors, data_plots$name),
    "trees" = setNames(e$listTree, data_plots$name),
    "species" = setNames(e$listSp, data_plots$name)
  )
  
  return(db_list)
}
