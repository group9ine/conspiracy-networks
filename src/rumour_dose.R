source("src/utils.R")

# rumour spreading as in Pastor-Satorass' review, but with doses
rumour_dose <- function(
  graph, n_iters, inf_0,
  p_skep, spr_rate, rec_rate, thresh,
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

  doses <- rep(0, n_nodes)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)
  k_inf <- rep(0, n_iters)

  when_inf <- rep(-1, n_nodes)
  reached <- rep(FALSE, n_nodes)
  dir_rec <- rep(TRUE, n_nodes)

  pick <- if (inf_0) inf_0 else floor(runif(1, 1, n_nodes + 1))
  inf[pick] <- TRUE
  reached[pick] <- TRUE
  doses[pick] <- thresh
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
    new_inf <- !rec & over
    inf <- (inf & !under) | new_inf
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
    k_inf[t] <- mean(k[new_inf])

    if (save_plots) {
      cols[inf] <- sir_pal[4]
      cols[rec] <- sir_pal[5]
      jpeg(sprintf("img/vid/dose_%04i.jpg", t), 1000, 1000, bg = sir_pal[1])
      plot(
        graph, layout = graph_lay,
        vertex.color = cols, vertex.size = vsize,
        vertex.label = NA, vertex.frame.color = NA,
        edge.color = sir_pal[2], edge.width = ewidth
      )
      dev.off()
    }

    if (n_inf[t] == 0) {
      n_sus <- n_sus[1:t]
      n_inf <- n_inf[1:t]
      n_rec <- n_rec[1:t]
      k_inf <- k_inf[1:t]
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
    doses = doses, sus = sus, rec = rec, k_inf = k_inf,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
