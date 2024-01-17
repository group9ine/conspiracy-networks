library(igraph)
library(ggplot2)
source("src/utils.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

dl2 <- read.csv("data/fb-orgs/L2.csv", col.names = c("from", "to"))
gl2 <- graph_from_data_frame(dl2, directed = FALSE) |> simplify()
V(gl2)$name <- seq_len(vcount(gl2))  # normalize naming

k <- unname(degree(gl2))
plot(logbins(k))
hist(k)
k

contagion <- function(graph, n_iters, n_inf, c_rate, d_wind, thresh) {
  degs <- igraph::degree(graph)
  n_nodes <- igraph::vcount(graph)
  ##[1]
  tol <- 5e-3 * n_nodes

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
      nbs <- nbs[degs[i] * runif(length(nbs)) < c_rate] ## replace with [2]
      doses[nbs, d_wind] <- doses[nbs, d_wind] + 1
    }

    inf <- rowSums(doses) >= thresh
    prev[t] <- sum(inf)

    if (t > 30 && sd(prev[seq(t - 30, t - 1)]) < tol)
      break
  }

  return(prev / n_nodes)
}

##proposed changes
## c_rate extracted from certain pdf, capped at c_rate=degs[i]
##[1] c_rate <- rpois(n=n_nodes,lambda=c_rate_lambda)
##[2] degs[i] * runif(length(nbs)) < min(degs[i],c_rate[i])
## optional plot function inside the main function
## plot = FALSE
##[3] if(plot == TRUE) {plot(prev / n_nodes)}
