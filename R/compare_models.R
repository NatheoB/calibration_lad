compare_models <- function(models_llpointwise) {
  
  # # Compute LOOC-CV and WAIC from the matrix of log-likelihood
  # waic_result <- tryCatch({
  #   loo::waic(pointwise_ll)
  # }, error = function(e) {
  #   message(e)
  #   NULL
  # })
  # 
  # loo_result <- tryCatch({
  #   loo::loo(pointwise_ll)
  # }, error = function(e) {
  #   message(e)
  #   NULL
  # })
  # 
  # 
  # return(list(
  #   "WAIC" = waic_result,
  #   "LOO-CV" = loo_result
  # ))
  
}