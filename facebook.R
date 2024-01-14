library(igraph)
library(ggplot2)
source("src/utils.R")

g <- read.csv("data/facebook-wall/graph.csv") |>
  graph_from_data_frame()
summary(g)
vcount(g)
ecount(g)

k <- unname(degree(g))

fit <- logbins(k)
data <- data.frame(x = seq_len(max(k)), y = tabulate(k))
ggplot(log10(data[data$y > 0, ]), aes(x, y)) +
  geom_point(alpha = 0.1) +
  geom_line(data = fit)
