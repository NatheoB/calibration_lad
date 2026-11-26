diagnose_convergence <- function(models_output_list, n_burning_per_rep) {
  
  # Convergence diagnostics of each model on the full post-burn-in chains (rank-normalised split-Rhat, bulk and tail ESS; Vehtari et al. 2021)
  # DREAMzs internal chains of a run are not independent (shared archive), so two Rhat are given:
  #   - rhat_within_run: max over runs of the Rhat of the internal chains of each run (has each run mixed?)
  #   - rhat_all: Rhat over all the chains of all the runs (do independent runs agree? the real convergence test)
  models_output_list <- purrr::discard(models_output_list, ~isTRUE(.x$failed))
  ids_model <- purrr::map_dbl(models_output_list, ~.x$id_model)
  
  models_output_list %>% 
    split(ids_model) %>% 
    purrr::map(function(runs) {
      
      par_names <- runs[[1]]$par_names
      
      # Post-burn-in internal chains of each run, as posterior draws arrays (iterations x chains x parameters)
      draws_runs <- purrr::map(runs, function(run) {
        n_burning_per_chain <- ceiling(n_burning_per_rep / length(run$chains))
        chains <- purrr::map(run$chains, ~as.matrix(.x)[-seq_len(n_burning_per_chain), par_names, drop = FALSE])
        posterior::as_draws_array(simplify2array(chains) %>% aperm(c(1, 3, 2)))
      })
      
      # Rhat within each run (max over runs)
      rhat_within_run <- draws_runs %>% 
        purrr::map(~posterior::summarise_draws(.x, "rhat")$rhat) %>% 
        purrr::reduce(pmax)
      
      # Rhat and ESS over all the chains of all the runs
      draws_all <- purrr::reduce(draws_runs, posterior::bind_draws, along = "chain")
      diagnostics <- posterior::summarise_draws(draws_all, "rhat", "ess_bulk", "ess_tail") %>% 
        dplyr::rename(parameter = variable, rhat_all = rhat) %>% 
        dplyr::mutate(rhat_within_run = rhat_within_run, .after = parameter)
      
      list(
        "id_model" = runs[[1]]$id_model,
        "n_reps" = length(runs),
        "n_chains" = posterior::nchains(draws_all),
        "n_iterations_per_chain" = posterior::niterations(draws_all),
        "diagnostics" = diagnostics,
        "summary" = diagnostics %>% dplyr::summarise(
          max_rhat_within_run = max(rhat_within_run),
          max_rhat_all = max(rhat_all),
          n_rhat_all_above_1.01 = sum(rhat_all > 1.01),
          min_ess_bulk = min(ess_bulk),
          min_ess_tail = min(ess_tail)
        )
      )
    })
}
