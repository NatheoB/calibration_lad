compare_models <- function(models_summary_pointwise) {
  
  # Pool the pointwise log-likelihood of all the runs of each model (failed runs removed) ----
  models_summary_pointwise <- purrr::discard(models_summary_pointwise, ~isTRUE(.x$failed))
  ids_model <- purrr::map_dbl(models_summary_pointwise, ~.x$id_model)
  
  models_loo <- models_summary_pointwise %>% 
    split(ids_model) %>% 
    purrr::map(function(x) {
      
      ll <- do.call(rbind, purrr::map(x, ~.x$llpointwise))
      sensors <- x[[1]]$sensors
      
      # PSIS-LOO (main criterion) and WAIC (secondary)
      # r_eff not given: draws thinned from long chains, assumed close to independent
      loo_result <- tryCatch(loo::loo(ll), error = function(e) {message(e); NULL})
      waic_result <- tryCatch(loo::waic(ll), error = function(e) {message(e); NULL})
      
      # Grouped PSIS-LOO by site (diagnostic of site-level dependence, many high Pareto k expected)
      ll_site <- t(rowsum(t(ll), group = sensors$site))
      loo_site <- tryCatch(loo::loo(ll_site), error = function(e) {message(e); NULL})
      
      # Sensors with an unreliable PSIS-LOO estimate (Pareto k > 0.7)
      sensors_highk <- NULL
      if (!is.null(loo_result)) {
        sensors_highk <- sensors %>% 
          dplyr::mutate(pareto_k = loo::pareto_k_values(loo_result)) %>% 
          dplyr::filter(pareto_k > 0.7)
      }
      
      list(
        "id_model" = x[[1]]$id_model,
        "n_reps" = length(x),
        "n_draws" = nrow(ll),
        "loo" = loo_result,
        "waic" = waic_result,
        "loo_site" = loo_site,
        "sensors_highk" = sensors_highk
      )
    })
  
  # Paired comparison of the models (elpd_diff and se_diff relative to the best model) ----
  loos <- purrr::compact(purrr::map(models_loo, "loo"))
  comparison <- NULL
  if (length(loos) > 1) comparison <- loo::loo_compare(loos)
  
  return(list(
    "models" = models_loo,
    "comparison" = comparison
  ))
}
