source("src/utils.R")

# rumour spreading without doses
rumour_base <- function(
  graph, n_iters, inf_0,
  p_skep, spr_rate, rec_rate,
  seed = FALSE, display = FALSE
) {
  if (seed) set.seed(seed)

  # we need to access by row later, so convert to row-ordered sp. mat.
  adj_mat <- as(igraph::as_adj(graph), "RsparseMatrix")
  n_nodes <- igraph::vcount(graph)

  sus <- rep(TRUE, n_nodes)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_sus <- rep(0, n_iters)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  reached <- rep(FALSE, n_nodes)
  dir_rec <- rep(FALSE, n_nodes)

  # pick a node randomly and infect it
  pick <- sample.int(n_nodes, inf_0, replace = FALSE)
  sus[pick] <- FALSE
  inf[pick] <- TRUE
  reached[pick] <- TRUE

  tol <- 1e-3 * n_nodes  # for the stopping condition
  for (t in seq_len(n_iters)) {
    # loop over infected nodes
    for (i in seq_len(n_nodes)[inf]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      s_mask <- sus[nbs]  # mask for every susceptible neighbour
      # subset the sus neighbors with probability spr_rate
      s_nbs <- nbs[s_mask][runif(sum(s_mask)) < spr_rate]
      # on this subset, select nodes that will recover with prob. p_skep
      r_mask <- runif(length(s_nbs)) < p_skep

      # i can recover with probability rec_rate for each contact
      # with an infected or recovered neighbour
      is_rec <- any(runif(sum(inf[nbs]) + sum(rec[nbs])) < rec_rate)
      inf[i] <- !is_rec
      rec[i] <- is_rec

      # update neighbours
      sus[s_nbs] <- FALSE
      rec[s_nbs[r_mask]] <- TRUE
      inf[s_nbs[!r_mask]] <- TRUE

      dir_rec[s_nbs[r_mask]] <- TRUE
    }


    # update counter vectors
    n_sus[t] <- sum(sus)
    n_inf[t] <- sum(inf)
    n_rec[t] <- sum(rec)
    reached <- reached | inf

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
    start = pick, reached = reached, dir_rec = dir_rec,
    sus = sus, rec = rec, inf = inf,
    n_sus = n_sus, n_inf = n_inf, n_rec = n_rec
  ))
}
