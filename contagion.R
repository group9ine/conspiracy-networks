library(igraph)
library(ggplot2)
source("src/utils.R")

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE)

k <- unname(degree(gl2))
data <- k

contagion <- function(graph, n_iters, post_prob, dose_window, thresh) {
  adj_mat <- igraph::get.adjacency(graph)
  n_nodes <- igraph::vcount(graph)

  doses <- matrix(0, nrow = n_nodes, ncol = dose_window)
  inf <- rep(FALSE, n_nodes)
  prev <- rep(0, n_iters)

  for (t in seq_len(dose_window)) {
    for (i in seq_len(n_nodes)[inf]) {
      if (runif(1) < post_prob) {
        nbs <- neighbors(graph, i, mode = "out")
        doses[nbs, t] <- doses[nbs, t] + 1
      }
    }
  }

  inf <- rowSums(doses) > thresh
  prev[dose_window] <- sum(inf)

  for (t in seq(dose_window, n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    for (i in seq_len(n_nodes)[inf]) {
      if (runif(1) < post_prob) {
        nbs <- neighbors(graph, i, mode = "out")
        doses[nbs, dose_window] <- doses[nbs, dose_window] + 1
      }
    }

    inf <- rowSums(doses) > thresh
    prev[t] <- sum(inf)
  }

  return(prev)
}
