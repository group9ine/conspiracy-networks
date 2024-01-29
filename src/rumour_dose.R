source("src/utils.R")

# rumour spreading as in Pastor-Satorass' review, but with doses
rumour_dose <- function(
  graph, n_iters, inf_0,
  p_skep, spr_rate, rec_rate, thresh,
  seed = FALSE, display = FALSE
) {
  if (seed) set.seed(seed)

  n_nodes <- igraph::vcount(graph)
  adj_mat <- as(
    igraph::as_adj(
      graph, attr = if (igraph::is_weighted(graph)) "weight" else NULL
    ),
    "RsparseMatrix"
  )

  doses <- rep(0, n_nodes)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  when_inf <- rep(-1, n_nodes)
  reached <- rep(FALSE, n_nodes)
  dir_rec <- rep(TRUE, n_nodes)

  pick <- if (inf_0) inf_0 else sample.int(n_nodes, 1, replace = FALSE)
  inf[pick] <- TRUE
  reached[pick] <- TRUE
  doses[pick] <- thresh[pick]
  when_inf[pick] <- 0

  sus <- !(inf | rec)

  for (t in seq_len(n_iters)) {
    # loop over infected nodes
    for (i in seq_len(n_nodes)[inf]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      crt <- adj_mat[i,, drop = FALSE]@x
      s_mask <- sus[nbs]  # mask for every sus neighbour

      # i can lose dose with probability rec_rate * c_rate for each contact
      # with an infected or recovered neighbour
      is_rec <- any(runif(sum(!s_mask)) < rec_rate * crt[!s_mask])
      doses[i] <- doses[i] - 1 * is_rec

      # subset the sus neighbours with probability spr_rate * c_rate
      c_mask <- runif(sum(s_mask)) < spr_rate * crt[s_mask]
      s_nbs <- nbs[s_mask][c_mask]
      # on this subset, gain or lose dose with probability p_skep
      doses[s_nbs] <- doses[s_nbs] + 1 - 2 * (runif(sum(c_mask)) < p_skep)
    }

    # update boolean vectors
    over <- doses >= thresh
    under <- doses < 0
    inf <- (inf & !under) | (!rec & over)
    rec <- rec | under
    sus <- !(inf | rec)
    when_inf[!reached & inf] <- t
    reached <- reached | inf
    # dir_rec starts from all T, goes F where a node goes above the thr.
    dir_rec <- dir_rec & !over

    # update counter vectors
    n_sus[t] <- sum(sus)
    n_inf[t] <- sum(inf)
    n_rec[t] <- sum(rec)

    if (n_inf[t] == 0) {
      n_sus <- n_sus[1:t]
      n_inf <- n_inf[1:t]
      n_rec <- n_rec[1:t]
      break
    }
  }

  # subset dir_rec with the actually recovered (leaving out those
  # who are still sus or inf)
  dir_rec <- dir_rec & rec
  n_sus <- n_sus / n_nodes
  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes

  if (display) print(ggsir(n_inf, n_rec, n_sus))

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(list(
    start = pick, duration = t, when_inf = when_inf,
    reached = reached, dir_rec = dir_rec,
    doses = doses, sus = sus, rec = rec,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
