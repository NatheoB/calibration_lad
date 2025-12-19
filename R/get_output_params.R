get_output_params <- function(models_output_list,
                              n_analysis,
                              thinning) {
  
  models_output_list %>%
    
    # Group list elements by models id 
    split(map_chr(., ~as.character(.x$id_model))) %>% 
    
    # Get the given sample 
    purrr::map(~purrr::map_dfr(.x, function(outs) {
      outs$outputs %>%
        getSample(start = nrow(outs$outputs$chain[[1]]) - round(n_analysis/length(outs$outputs$chain)),
                  end = NULL,
                  thin = thinning,
                  coda = TRUE) %>%
        purrr::map(function(out) {
          as.data.frame(out) %>%
            dplyr::mutate(iteration = row_number()) %>%
            dplyr::relocate(iteration)
        }) %>%
        dplyr::bind_rows(.id = "subchain") %>% 
        dplyr::mutate(chain = outs$i_chain, .before = subchain) %>% 
        dplyr::mutate(id_set = paste(chain, subchain, iteration, sep = "_"),
                      .before = chain)
    }))
  
}