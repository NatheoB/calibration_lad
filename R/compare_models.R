compare_models <- function(models_summary_pointwise) {
  
  models_summary_pointwise %>% 
    purrr::map(function(x) {
      
      # Compute LOOC-CV and WAIC from the matrix of log-likelihood
      waic_result <- tryCatch({
        loo::waic(x$llpointwise)$estimates
      }, error = function(e) {
        message(e)
        NULL
      })
      
      loo_result <- tryCatch({
        loo::loo(x$llpointwise)$estimates
      }, error = function(e) {
        message(e)
        NULL
      })
      
      return(list(
        "WAIC" = waic_result,
        "LOOCV" = loo_result
      ))
      
    })
  
}