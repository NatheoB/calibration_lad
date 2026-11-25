evaluate_models <- function(models_summary_pointwise) {
  
  data_residuals <- models_summary_pointwise %>%
    purrr::map(
      ~.x$residuals %>% 
        as.data.frame() %>% 
        dplyr::mutate(iteration = row_number()) %>%
        tidyr::pivot_longer(!iteration,
                            names_to = "id_sensor",
                            values_to = "residuals",
                            names_prefix = "V",
                            names_transform = as.numeric) %>% 
        dplyr::mutate(
          id_model = .x$id_model,
          i_rep = .x$i_rep,
          .before = iteration
        )
    )
  
  evaluation_indicators <- data_residuals %>% 
    purrr::map(
       ~.x %>% 
        dplyr::group_by(id_model, i_rep, iteration) %>% 
        dplyr::summarise(
          MAE = mean(abs(residuals)),
          RMSE = sqrt(mean(residuals^2))
        )
    )
  
  # Predictions per sensor and model (runs pooled) ----
  models_summary_pointwise <- purrr::discard(models_summary_pointwise, ~isTRUE(.x$failed))
  ids_model <- purrr::map_dbl(models_summary_pointwise, ~.x$id_model)
  
  sensors_pred <- models_summary_pointwise %>% 
    split(ids_model) %>% 
    purrr::map(function(x) {
      
      res <- do.call(rbind, purrr::map(x, ~.x$residuals))
      sigma <- unlist(purrr::map(x, ~.x$sigma))
      
      # Posterior predictive draws minus observed PACL (residual + observation noise of each draw)
      res_pred <- res + matrix(rnorm(length(res), 0, sigma), nrow = nrow(res))
      
      x[[1]]$sensors %>% 
        dplyr::mutate(
          id_model = x[[1]]$id_model,
          pacl_pred_mean = pacl_field + colMeans(res),
          pacl_pred_q05 = pacl_field + apply(res_pred, 2, quantile, 0.05),
          pacl_pred_q25 = pacl_field + apply(res_pred, 2, quantile, 0.25),
          pacl_pred_q75 = pacl_field + apply(res_pred, 2, quantile, 0.75),
          pacl_pred_q95 = pacl_field + apply(res_pred, 2, quantile, 0.95)
        )
    }) %>% 
    dplyr::bind_rows()
  
  # Summary per model: error of the posterior mean prediction and coverage of predictive intervals ----
  summary_models <- sensors_pred %>% 
    dplyr::group_by(id_model) %>% 
    dplyr::summarise(
      MAE = mean(abs(pacl_pred_mean - pacl_field)),
      RMSE = sqrt(mean((pacl_pred_mean - pacl_field)^2)),
      bias = mean(pacl_pred_mean - pacl_field),
      coverage_50 = mean(pacl_field >= pacl_pred_q25 & pacl_field <= pacl_pred_q75),
      coverage_90 = mean(pacl_field >= pacl_pred_q05 & pacl_field <= pacl_pred_q95)
    )
  
  return(list(
    "residuals" = data_residuals,
    "indicators" = evaluation_indicators,
    "sensors_pred" = sensors_pred,
    "summary_models" = summary_models
  ))
}