save_output_models <- function(exp_design, 
                               models_setup, 
                               models_output,
                               models_summary_pointwise,
                               models_comparison,
                               models_evaluation,
                               output_params,
                               output_folder,
                               output_name) {
  
  output_models_list <- list(
    "exp_design" = exp_design,
    "setup" = models_setup,
    "output" = models_output,
    "pointwise" = models_summary_pointwise,
    "comparison" = models_comparison,
    "evaluation" = models_evaluation,
    "params" = output_params
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_models_list, file = file.path(output_folder, output_name))
}