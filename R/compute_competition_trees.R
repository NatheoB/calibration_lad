compute_competition_trees <- function(data_trees, 
                                      local_ba_radius,
                                      plot_area, 
                                      plot_xsize, 
                                      plot_ysize) {
  
  # ── Plot weight (for per-hectare scaling) ──────────────────────────────────
  plot_weight <- 1 / plot_area
  local_weight <- 10000 / (pi * local_ba_radius^2)
  
  # ── Stand-level competition variables ─────────────────────────────────────
  data_trees <- data_trees %>% 
    dplyr::mutate(
      dg_cm         = sqrt(sum(dbh_cm^2) / n()),
      compet_status = dbh_cm / dg_cm,
      G_m2          = (pi * (dbh_cm / 100)^2 / 4),
      G_m2ha        = G_m2 * plot_weight,
      batot_m2ha    = sum(G_m2ha)
    ) %>%
    dplyr::arrange(dplyr::desc(dbh_cm)) %>%
    dplyr::mutate(bal_m2ha = ave(G_m2ha, FUN = cumsum) - G_m2ha) %>%
    dplyr::arrange(dplyr::desc(h_m)) %>%
    dplyr::mutate(bah_m2ha = ave(G_m2ha, FUN = cumsum) - G_m2ha)
  
  # ── Toroidal distance helper ───────────────────────────────────────────────
  # The shortest distance
  torus_dist <- function(xi, xj, size) {
    d <- xj - xi
    d - size * round(d / size)   # equivalent to modular wrap
  }
  
  # ── Local indices (computed row-by-row via sapply) ─────────────────────────
  n_trees <- nrow(data_trees)
  
  batot_local   <- numeric(n_trees)
  bal_local     <- numeric(n_trees)
  bah_local     <- numeric(n_trees)
  
  for (i in seq_len(n_trees)) {
    
    dx <- torus_dist(data_trees$x[i], data_trees$x, plot_xsize)
    dy <- torus_dist(data_trees$y[i], data_trees$y, plot_ysize)
    dist_to_i <- sqrt(dx^2 + dy^2)
    
    # Neighbours: within radius, excluding focal tree itself
    neighbours <- which(dist_to_i <= local_ba_radius & dist_to_i > 0)
    
    if (length(neighbours) == 0) {
      batot_local[i]  <- 0
      bal_local[i]    <- 0
      bah_local[i]    <- 0
    } else {
      G_neighbours     <- data_trees$G_m2[neighbours] * local_weight
      batot_local[i]   <- sum(G_neighbours)
      bal_local[i]     <- sum(G_neighbours[data_trees$dbh_cm[neighbours] > data_trees$dbh_cm[i]])
      bah_local[i]     <- sum(G_neighbours[data_trees$h_m[neighbours]   > data_trees$h_m[i]])
    }
  }
  
  # ── Attach local indices and return ───────────────────────────────────────
  data_trees %>%
    dplyr::mutate(
      batot_local_m2ha  = batot_local,
      bal_local_m2ha    = bal_local,
      bah_local_m2ha    = bah_local
    )
}
