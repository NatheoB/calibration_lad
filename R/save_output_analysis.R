save_output_analysis <- function(output_params,
                                 output_tree_lad,
                                 output_tree_light,
                                 output_folder,
                                 output_name) {
  
  output_analysis_list <- list(
    "output_params" = output_params,
    "output_tree_lad" = output_tree_lad,
    "output_tree_light" = output_tree_light
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_analysis_list, file = file.path(output_folder, output_name))
}