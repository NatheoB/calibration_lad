get_output_params <- function(models_output_list,
                              n_analysis) {
  
  tmp_out <- models_output_list %>%
    
    # Group list elements by models id 
    split(purrr::map_chr(., ~as.character(.x$id_model))) %>% 
    
    # Get the given sample for a replicate and a model
    purrr::map(~purrr::map_dfr(.x, function(outs) {
      
      outs$outputs %>%
        
        BayesianTools::getSample(coda = FALSE,
                                 numSamples = n_analysis) %>% 
        as.data.frame() %>% 
        dplyr::mutate(rep = outs$i_rep, 
                      id_params = row_number(),
                      id_set = paste(rep, id_params, sep = "_")) %>% 
        dplyr::relocate(id_set, rep, id_params)
    }))
  
}