save_models <- function(exp_design, 
                        models_setup, 
                        models_output,
                        models_summary_pointwise,
                        models_comparison,
                        models_evaluation,
                        output_folder,
                        output_name) {
  
  out_list <- list(
    "exp_design" = exp_design,
    "setup" = models_setup,
    "output" = models_output,
    "pointwise" = models_summary_pointwise,
    "comparison" = models_comparison,
    "evaluation" = models_evaluation
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(out_list, file = file.path(output_folder, output_name))
}