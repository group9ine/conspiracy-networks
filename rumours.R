library(igraph)
library(ggplot2)
library(data.table)
source("src/utils.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE) |> simplify()
V(gl2)$name <- seq_len(vcount(gl2))  # normalize naming

##############
# COPENHAGEN #
##############

# phone calls
dc <- fread(
  "data/copenhagen/calls.csv",
  col.names = c("time", "from", "to", "duration")
)[duration > 0, .N, keyby = .(a = pmin(from, to), b = pmax(from, to))]
# sms
ds <- fread(
  "data/copenhagen/sms.csv",
  col.names = c("time", "from", "to")
)[, .N, keyby = .(a = pmin(from, to), b = pmax(from, to))]
# facebook friends
df <- fread("data/copenhagen/fb_friends.csv", col.names = c("a", "b"))

# big merge with custom weighting
dcph <- merge(dc, merge(ds, df, all = TRUE), by = c("a", "b"), all = TRUE)[
  is.na(N.x), N.x := 0][
    is.na(N.y), N.y := 0][
      N.x + N.y > 0, weight := 1 + log(N.x + N.y)][
        N.x + N.y == 0, weight := 1][
          , c("N.x", "N.y") := NULL]

gcph <- graph_from_data_frame(dcph, directed = FALSE) |> simplify()
V(gcph)$name <- seq_len(vcount(gcph))  # normalize naming
is_simple(gcph)
clusters(gcph)  # only one cc, good

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

  if (display) {
    plt <- data.frame(
      iter = seq_along(n_inf),
      inf = n_inf, rec = n_rec, sus = n_sus
    ) |>
      tidyr::pivot_longer(-iter, names_to = "class", values_to = "pop") |>
      dplyr::mutate(class = factor(class, levels = c("sus", "inf", "rec"))) |>
      ggplot(ggplot2::aes(iter, pop, colour = class)) +
        ggplot2::geom_line() +
        ggplot2::scale_colour_manual(
          values = c("darkgoldenrod1", "firebrick", "dodgerblue4"),
          labels = c("Susceptible", "Infected", "Recovered")
        ) +
        ggplot2::labs(
          x = "Time step",
          y = "Fractional prevalence",
          colour = "Compartment"
        )
    print(plt)
  }

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(data.frame(sus = n_sus, inf = n_inf, rec = n_rec))
}

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
    # loop over infected nodes
    for (i in seq_len(n_nodes)[inf]) {
      # this returns the indices of the i-th row's non-zero elements
      # i.e. out-neighbors of node i
      nbs <- adj_mat[i,, drop = FALSE]@j + 1
      nmask <- runif(length(nbs)) < p_talk
      nbs <- nbs[nmask]
      wgt <- adj_mat[i,, drop = FALSE]@x[nmask]

      smask <- sus[nbs]
      snbs <- nbs[smask]

      # susceptible neighbors lose or gain doses with prob. = p_skep
      doses[snbs, d_wind] <- doses[snbs, d_wind] +
        (1 - 2 * (runif(sum(smask)) < p_skep)) * wgt[smask]
      # current node (infected) loses doses for each neighbouring I or R
      doses[i, d_wind] <- doses[i, d_wind] - sum(runif(sum(!smask)) < p_stop)
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

  if (display) {
    plt <- data.frame(
      iter = seq_along(n_inf),
      inf = n_inf, rec = n_rec, sus = n_sus
    ) |>
      tidyr::pivot_longer(-iter, names_to = "class", values_to = "pop") |>
      dplyr::mutate(class = factor(class, levels = c("sus", "inf", "rec"))) |>
      ggplot(ggplot2::aes(iter, pop, colour = class)) +
        ggplot2::geom_line() +
        ggplot2::scale_colour_manual(
          values = c("darkgoldenrod1", "firebrick", "dodgerblue4"),
          labels = c("Susceptible", "Infected", "Recovered")
        ) +
        ggplot2::labs(
          x = "Time step",
          y = "Fractional prevalence",
          colour = "Compartment"
        )
    print(plt)
  }

  # clean up the rng seed if it was set
  rm(.Random.seed, envir = .GlobalEnv)

  return(data.frame(sus = n_sus, inf = n_inf, rec = n_rec))
}

# null model
ger <- sample_gnm(n = vcount(gcph), m = ecount(gcph))
res_er <- rumour_skep(
  graph = ger, n_iters = 100, inf_0 = 1,
  p_talk = 0.8, p_skep = 0.2, p_stop = 0.1,
  d_wind = 7, thresh = 5,
  seed = FALSE, display = TRUE
)
res_cph <- rumour_skep(
  graph = gcph, n_iters = 100, inf_0 = 1,
  p_talk = 0.4, p_skep = 0.2, p_stop = 0.1,
  d_wind = 7, thresh = 5,
  seed = FALSE, display = TRUE
)
