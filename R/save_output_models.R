save_output_models <- function(exp_design, 
                               models_setup, 
                               models_output,
                               output_folder,
                               output_name) {
  
  output_models_list <- list(
    "exp_design" = exp_design,
    "setup" = models_setup,
    "output" = models_output
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_models_list, file = file.path(output_folder, output_name))
}