source("src/utils.R")

# rumour spreading with random contact rates (log-normally distributed)
rumour_rnd <- function(
  graph, n_iters, inf_0, c_rate_mu, c_rate_sig, d_wind, thresh,
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
  # randomly assign a contact rate to each node
  rates <- rlnorm(n_nodes, c_rate_mu, c_rate_sig)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  pick <- sample.int(n_nodes, inf_0, replace = FALSE)
  inf[pick] <- TRUE
  doses[pick, d_wind] <- thresh

  tol <- 1e-3 * n_nodes  # for the stopping condition
  for (t in seq_len(n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    # loop over infected nodes
    for (i in seq_len(n_nodes)[inf]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      wgt <- adj_mat[i,, drop = FALSE]@x
      mask <- degs[i] * runif(length(nbs)) < rates[i]
      nbs <- nbs[mask]
      doses[nbs, d_wind] <- doses[nbs, d_wind] + wgt[mask]
    }

    # update boolean vectors
    over <- rowSums(doses) >= thresh
    new_inf <- !rec & over
    rec <- rec | (inf & !over)
    inf <- new_inf

    # update counter vectors
    n_inf[t] <- sum(inf)
    n_rec[t] <- sum(rec)

    if (t > 30 && sd(n_inf[(t - 30):t]) < tol) {
      n_inf <- n_inf[1:t]
      n_rec <- n_rec[1:t]
      break
    }
  }

  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes
  n_sus <- 1 - n_inf - n_rec

  if (display) print(ggsir(n_inf, n_rec, n_sus))

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(data.frame(sus = n_sus, inf = n_inf, rec = n_rec))
}
