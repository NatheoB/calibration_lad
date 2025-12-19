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
          i_chain = .x$i_chain,
          .before = iteration
        )
    )
  
  evaluation_indicators <- data_residuals %>% 
    purrr::map(
       ~.x %>% 
        dplyr::group_by(id_model, i_chain, iteration) %>% 
        dplyr::summarise(
          MAE = mean(abs(residuals)),
          RMSE = sqrt(mean(residuals^2))
        )
    )
  
  return(list(
    "residuals" = data_residuals,
    "indicators" = evaluation_indicators
  ))
}