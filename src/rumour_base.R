source("src/utils.R")

# rumour spreading without doses
rumour_base <- function(
  graph, n_iters, inf_0,
  p_skep, spr_rate, rec_rate,
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

  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  when_inf <- rep(-1, n_nodes)
  reached <- rep(FALSE, n_nodes)
  dir_rec <- rep(FALSE, n_nodes)

  # pick a node randomly and infect it
  pick <- if (inf_0) inf_0 else sample.int(n_nodes, 1, replace = FALSE)
  inf[pick] <- TRUE
  reached[pick] <- TRUE
  when_inf[pick] <- 0

  sus <- !(inf | rec)

  for (t in seq_len(n_iters)) {
    # loop over infected nodes
    new_inf <- rep(FALSE, n_nodes)
    new_rec <- rep(FALSE, n_nodes)
    for (i in seq_len(n_nodes)[inf]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      crt <- adj_mat[i,, drop = FALSE]@x
      s_mask <- sus[nbs]  # mask for every sus neighbour

      # i can recover with probability rec_rate * c_rate for each contact
      # with an infected or recovered neighbour
      is_rec <- any(runif(sum(!s_mask)) < rec_rate * crt[!s_mask])
      new_inf[i] <- !is_rec
      new_rec[i] <- is_rec

      # subset the sus neighbors with probability spr_rate * c_rate
      s_nbs <- nbs[s_mask][runif(sum(s_mask)) < spr_rate * crt[s_mask]]
      # on this subset, select nodes that will recover with prob. p_skep
      r_mask <- runif(length(s_nbs)) < p_skep

      # update neighbours
      new_rec[s_nbs[r_mask]] <- TRUE
      new_inf[s_nbs[!r_mask]] <- TRUE

      dir_rec[s_nbs[r_mask]] <- TRUE
    }

    inf <- (inf & !new_rec) | new_inf
    rec <- rec | new_rec
    sus <- !(inf | rec)
    when_inf[!reached & inf] <- t
    reached <- reached | inf

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

  n_sus <- n_sus / n_nodes
  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes

  if (display) print(ggsir(n_inf, n_rec, n_sus))

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(list(
    start = pick, duration = t,
    when_inf = when_inf, reached = reached,
    dir_rec = dir_rec, sus = sus, rec = rec,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
