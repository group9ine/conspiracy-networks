library(tidyverse)
library(igraph)

g <- read_csv("data/facebook-wall/graph.csv", show_col_types = FALSE) |>
  graph_from_data_frame()

head(E(g)$time)  # time is read correctly

# we can also use actual date + times:
g <- read_csv("data/facebook-wall/graph.csv", show_col_types = FALSE) |>
  mutate(time = as_datetime(time, tz = "UTC")) |>
  graph_from_data_frame()
head(E(g)$time)

times <- E(g)$time

m <- interval(first(times), last(times)) |>
  time_length("months")
sprintf("%d years and %d months", floor(m / 12), floor((m - floor(m)) * 12))
