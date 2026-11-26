get_output_params <- function(models_output_list,
                              n_burning,
                              n_samples_per_chain) {
  
  models_output_list %>%
    
    # Group list elements by models id 
    split(purrr::map_chr(., ~as.character(.x$id_model))) %>% 
    
    # Get the given sample for a replicate and a model
    purrr::map(~purrr::map_dfr(.x, function(outs) {
      
      # Samples of the internal chains of the run, after burn-in
      chain_samples <- get_chain_samples(outs, n_burning, n_samples_per_chain)
      
      chain_samples$samples %>%
        as.data.frame() %>% 
        
        # Rename parameters
        dplyr::rename_all(~paste0("p_", .)) %>% 
        
        # Set unique id to paremeters
        dplyr::mutate(rep = outs$i_rep, 
                      id_chain = chain_samples$id_chain,
                      id_params = row_number(),
                      id_set = paste(rep, id_params, sep = "_")) %>% 
        dplyr::relocate(id_set, rep, id_chain, id_params)
    }))
  
}