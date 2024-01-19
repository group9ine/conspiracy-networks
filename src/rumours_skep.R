source("src/utils.R")

# rumour spreading as in Pastor-Satorass' review
rumour_skep <- function(
  graph, n_iters, inf_0, p_talk, p_skep, p_stop, d_wind, thresh,
  seed = FALSE, display = FALSE
) {
  if (seed) set.seed(seed)

  # we need to access by row later, so convert to row-ordered sp. mat.
  adj_mat <- as(
    igraph::as_adj(graph, attr = if (is_weighted(graph)) "weight" else NULL),
    "RsparseMatrix"
  )
  degs <- igraph::degree(graph)
  n_nodes <- igraph::vcount(graph)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  pick <- sample.int(n_nodes, inf_0, replace = FALSE)
  inf[pick] <- TRUE
  doses[pick, d_wind] <- thresh

  sus <- !(inf | rec)

  tol <- 1e-3 * n_nodes  # for the stopping condition
  for (t in seq_len(n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    # loop over nodes that are both infected and active (with prob. p_talk)
    active <- runif(n_nodes) < p_talk
    for (i in seq_len(n_nodes)[inf & active]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      a_mask <- active[nbs]
      a_nbs <- nbs[a_mask]
      wgt <- adj_mat[i,, drop = FALSE]@x[a_mask]

      s_mask <- sus[a_nbs]
      s_nbs <- a_nbs[s_mask]

      # susceptible neighbors lose or gain doses with prob. = p_skep
      doses[s_nbs, d_wind] <- doses[s_nbs, d_wind] +
        (1 - 2 * (runif(sum(s_mask)) < p_skep)) * wgt[s_mask]
      # current node (infected) loses doses for each neighbouring I or R
      doses[i, d_wind] <- doses[i, d_wind] - sum(runif(sum(!s_mask)) < p_stop)
    }

    # update boolean vectors
    over <- rowSums(doses) >= thresh
    new_inf <- !rec & over
    rec <- rec | (inf & !over)
    inf <- new_inf
    sus <- !(inf | rec)

    # update counter vectors
    n_sus[t] <- sum(sus)
    n_inf[t] <- sum(inf)
    n_rec[t] <- sum(rec)

    if (t > 30 && sd(n_inf[(t - 30):t]) < tol) {
      n_sus <- n_sus[1:t]
      n_inf <- n_inf[1:t]
      n_rec <- n_rec[1:t]
      break
    }
  }

  n_sus <- n_sus / n_nodes
  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes

  if (display) print(ggsir(n_inf, n_rec, n_sus))

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(list(
    sus = sus, rec = rec, inf = inf,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
