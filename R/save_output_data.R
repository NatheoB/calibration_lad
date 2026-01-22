save_output_data <- function(data_output_tree,
                             data_output_stand,
                             data_output_light,
                             output_folder,
                             output_name) {
  
  output_data_list <- list(
    "data_output_tree" = data_output_tree,
    "data_output_stand" = data_output_stand,
    "data_output_light" = data_output_light
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_data_list, file = file.path(output_folder, output_name))
}