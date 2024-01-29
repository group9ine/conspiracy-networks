source("src/utils.R")

# rumour spreading without doses
rumour_base <- function(
  graph, n_iters, inf_0,
  p_skep, spr_rate, rec_rate,
  seed = FALSE, display = FALSE,
  save_plots = FALSE, graph_lay = NULL
) {
  if (seed) set.seed(seed)

  n_nodes <- igraph::vcount(graph)
  adj_mat <- as(
    igraph::as_adj(
      graph, attr = if (igraph::is_weighted(graph)) "weight" else NULL
    ),
    "RsparseMatrix"
  )
  k <- igraph::degree(graph, mode = "out")

  if (save_plots) {
    cols <- rep(sir_pal[3], n_nodes)
    vsize <- 2 * k^0.3
    ewidth <- 0.5 * E(graph)$weight / min(E(graph)$weight)
  }

  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)
  k_inf <- rep(0, n_iters)

  when_inf <- rep(-1, n_nodes)
  reached <- rep(FALSE, n_nodes)
  dir_rec <- rep(FALSE, n_nodes)

  # pick a node randomly and infect it
  pick <- if (inf_0) inf_0 else floor(runif(1, 1, n_nodes + 1))
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
    k_inf[t] <- mean(k[new_inf])

    if (save_plots) {
      cols[inf] <- sir_pal[4]
      cols[rec] <- sir_pal[5]
      plot(
        graph, layout = graph_lay, edge.width = ewidth,
        vertex.color = cols, vertex.size = vsize,
        vertex.label = NA, vertex.frame.color = NA
      )

      file_name <- sprintf("plots/plot_%03i.png", t)
      dev.copy(png, file_name)
    }

    if (n_inf[t] == 0) {
      n_sus <- n_sus[1:t]
      n_inf <- n_inf[1:t]
      n_rec <- n_rec[1:t]
      k_inf <- k_inf[1:t]
      break
    }
  }

  n_sus <- n_sus / n_nodes
  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes

  if (display) print(ggsir(n_inf, n_rec, n_sus))
  if (save_plots) dev.off()

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(list(
    start = pick, duration = t,
    when_inf = when_inf, reached = reached,
    dir_rec = dir_rec, sus = sus, rec = rec, k_inf = k_inf,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
