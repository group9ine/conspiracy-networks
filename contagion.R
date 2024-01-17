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
  degs <- igraph::degree(graph)
  n_nodes <- igraph::vcount(graph)
  # randomly assign a contact rate to each node
  rates <- rlnorm(n_nodes, c_rate_mu, c_rate_sig)

  doses <- matrix(0, nrow = n_nodes, ncol = d_wind)
  # 1 = sus, 0 = inf, < 0 = rec
  status <- rep(1, n_nodes)
  n_inf <- rep(0, n_iters)
  n_rec <- rep(0, n_iters)

  pick <- sample.int(n_iters, inf_0, replace = FALSE)
  status[pick] <- 0
  doses[pick, d_wind] <- thresh

  tol <- 1e-3 * n_nodes  # for the stopping condition
  for (t in seq_len(n_iters)) {
    doses <- cbind(doses[, -1], rep(0, n_nodes))
    # loop over infected nodes
    for (i in seq_len(n_nodes)[!status]) {
      nbs <- igraph::neighbors(graph, i, mode = "out")
      nbs <- nbs[degs[i] * runif(length(nbs)) < rates[i]]
      doses[nbs, d_wind] <- doses[nbs, d_wind] + 1
    }

    status <- status - rowSums(doses) >= thresh
    n_inf[t] <- sum(!status)
    n_rec[t] <- sum(status < 0)

    #if (t > 30 && sd(n_rec[(t - 30):t]) < tol) {
    #  n_inf <- n_inf[1:t]
    #  n_rec <- n_rec[1:t]
    #  break
    #}
  }

  n_inf <- n_inf / n_nodes
  n_rec <- n_rec / n_nodes
  n_sus <- 1 - n_inf - n_rec

  if (display) {
    data.frame(
      iter = seq_along(n_inf), inf = n_inf, rec = n_rec, sus = n_sus
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
          x = "Time step", y = "Fractional prevalence", colour = "Compartment"
        )
  }

  return(data.frame(sus = n_sus, inf = n_inf, rec = n_rec))
}
