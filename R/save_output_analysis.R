save_output_analysis <- function(output_params,
                                 data_output_tree,
                                 data_output_stand,
                                 data_output_light,
                                 output_folder,
                                 output_name) {
  
  output_analysis_list <- list(
    "output_params" = output_params,
    "data_output_tree" = data_output_tree,
    "data_output_stand" = data_output_stand,
    "data_output_light" = data_output_light
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_analysis_list, file = file.path(output_folder, output_name))
}