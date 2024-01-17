library(igraph)
library(ggplot2)
source("src/utils.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE) |> simplify()
V(gl2)$name <- seq_len(vcount(gl2))  # normalize naming


contagion <- function(
  graph, n_iters, inf_0, c_rate_mu, c_rate_sig, d_wind, thresh, display = FALSE
) {
  # we need to access by row later, so convert to row-ordered sp. mat.
  adj_mat <- as(igraph::get.adjacency(graph), "RsparseMatrix")
  degs <- igraph::degree(graph)
  n_nodes <- igraph::vcount(graph)
  # randomly assign a contact rate to each node
  rates <- rlnorm(n_nodes, c_rate_mu, c_rate_sig)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  inf <- rep(FALSE, n_nodes)
  rec <- rep(FALSE, n_nodes)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  pick <- sample.int(n_iters, inf_0, replace = FALSE)
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
      nbs <- nbs[degs[i] * runif(length(nbs)) < rates[i]]
      doses[nbs, d_wind] <- doses[nbs, d_wind] + 1
    }

    # update boolean vectors
    over <- rowSums(doses) >= thresh
    new_inf <- !rec & over
    rec <- rec | (inf & !over)
    inf <- new_inf

    # update counter vectors
    n_inf[t] <- sum(inf)
    n_rec[t] <- sum(rec)
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

  return(data.frame(sus = n_sus, inf = n_inf, rec = n_rec))
}
