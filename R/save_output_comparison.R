save_output_comparison <- function(models_summary_pointwise,
                                   models_comparison,
                                   models_evaluation,
                                   output_folder,
                                   output_name) {
  
  output_comparison_list <- list(
    "pointwise" = models_summary_pointwise,
    "comparison" = models_comparison,
    "evaluation" = models_evaluation
  )
  
  dir.create(output_folder, showWarnings = FALSE)
  save(output_comparison_list, file = file.path(output_folder, output_name))
}