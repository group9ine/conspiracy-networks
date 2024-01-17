library(igraph)
library(ggplot2)
source("src/utils.R")

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE) |> simplify()
V(gl2)$name <- seq_len(vcount(gl2))  # normalize naming


contagion <- function(graph, n_iters, n_inf, c_rate, d_wind, thresh) {
  degs <- igraph::degree(graph)
  n_nodes <- igraph::vcount(graph)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  inf <- rep(FALSE, n_nodes)
  prev <- rep(0, n_iters)

  # infect n_inf nodes at random
  pick <- sample.int(n_iters, n_inf, replace = FALSE)
  inf[pick] <- TRUE
  doses[pick, d_wind] <- thresh

  for (t in seq_len(n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    for (i in seq_len(n_nodes)[inf]) {
      nbs <- igraph::neighbors(graph, i, mode = "out")
      nbs <- nbs[degs[i] * runif(length(nbs)) < c_rate]
      doses[nbs, d_wind] <- doses[nbs, d_wind] + 1
    }

    inf <- rowSums(doses) >= thresh
    prev[t] <- sum(inf)

    if (t > 30 && all(prev[seq(t - 30, t - 1)] == prev[t])) break
  }

  return(prev / n_nodes)
}
