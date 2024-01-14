library(igraph)
library(data.table)

fb <- fread(
  "data/facebook-wall/graph.tsv",
  sep = " ",
  col.names = c("from", "to", "count", "time")
)[rep(seq_len(.N), count), -"count"]  # duplicate rows based on 'count' column

g <- graph_from_data_frame(fb)
head(E(g)$time)  # time info is there, it seems
