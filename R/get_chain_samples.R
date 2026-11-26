get_chain_samples <- function(run_output, n_burning_per_rep, n_samples_per_rep) {
  
  # DREAMzs runs several internal chains in parallel (3 here) and BayesianTools counts iterations over all of them:
  # the burn-in and the number of samples of a run are thus split between the internal chains
  chains <- run_output$chains
  if (!coda::is.mcmc.list(chains)) stop("The run output must contain the internal chains as a coda::mcmc.list")
  n_chains <- length(chains)
  
  n_burning_per_chain <- ceiling(n_burning_per_rep / n_chains)
  n_samples_per_chain <- rep(ceiling(n_samples_per_rep / n_chains), n_chains) # Same number in every chain (needed by loo::relative_eff()): n_samples_per_rep is rounded up to a multiple of the number of chains
  
  # For each internal chain: remove the burn-in, then keep evenly spaced samples of the parameters
  samples_list <- lapply(seq_len(n_chains), function(i_chain) {
    
    chain <- as.matrix(chains[[i_chain]])[, run_output$par_names, drop = FALSE]
    if (nrow(chain) <= n_burning_per_chain) stop("Internal chain ", i_chain, " has ", nrow(chain), " iterations, not more than its burn-in (", n_burning_per_chain, ")")
    
    chain <- chain[-seq_len(n_burning_per_chain), , drop = FALSE]
    chain[round(seq(1, nrow(chain), length.out = n_samples_per_chain[i_chain])), , drop = FALSE]
  })
  
  list(
    "samples" = do.call(rbind, samples_list),
    "id_chain" = rep(seq_len(n_chains), times = n_samples_per_chain) # Internal chain of each sample (for r_eff and convergence diagnostics)
  )
}
