library(igraph)
library(ggplot2)
source("src/utils.R")

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE) |> simplify()
V(gl2)$name <- seq_len(vcount(gl2))  # normalize naming


contagion <- function(graph, n_iters, c_prob, d_wind, thresh) {
  adj_mat <- igraph::get.adjacency(graph)
  n_nodes <- igraph::vcount(graph)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  inf <- rep(FALSE, n_nodes)
  prev <- rep(0, n_iters)

  # infect a node at random
  pick <- floor(runif(1, min = 1, max = n_nodes + 1))
  inf[pick] <- TRUE
  doses[pick, d_wind] <- thresh

  for (t in seq(1, n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    for (i in seq_len(n_nodes)[inf]) {
      if (runif(1) < c_prob) {
        nbs <- adj_mat[, i, drop = FALSE]@i + 1
        doses[nbs, d_wind] <- doses[nbs, d_wind] + 1
      }
    }

    inf <- rowSums(doses) >= thresh
    prev[t] <- sum(inf)
  }

  return(prev / n_nodes)
}
