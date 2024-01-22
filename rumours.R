library(igraph)
library(ggplot2)
library(data.table)
source("src/utils.R")
source("src/rumour_base.R")
source("src/rumour_dose.R")

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
# dcph <- merge(dc, merge(ds, df, all = TRUE), all = TRUE)[
#   is.na(N.x), N.x := 0][is.na(N.y), N.y := 0][
#       N.x + N.y > 0, weight := 1 + log(N.x + N.y)][
#         N.x + N.y == 0, weight := 1][
#           , c("N.x", "N.y") := NULL]

# without facebook
dcph <- merge(dc, ds, all = TRUE)[
  is.na(N.x), N.x := 0][is.na(N.y), N.y := 0][
    , .(a, b, weight = 1 + log(N.x + N.y))]

gcph <- graph_from_data_frame(dcph, directed = FALSE) |>
  simplify() |>
  get_lcc()
V(gcph)$name <- seq_len(vcount(gcph))  # normalize naming
is_simple(gcph)
clusters(gcph)  # only one cc, good

k <- degree(gcph)
c <- closeness(gcph)
c_sc <- (c - min(c)) / (max(c) - min(c))
n_cols <- 20
V(gcph)$color <- viridis::plasma(n_cols + 1)[1 + floor(n_cols * c_sc)]
plot(
  gcph, layout = layout_with_graphopt(gcph),
  edge.width = 0.5 * E(gcph)$weight,
  vertex.size = 2 * k^0.3, vertex.label = NA, vertex.frame.color = NA
)

# null model
ger <- sample_gnm(n = vcount(gcph), m = ecount(gcph))
res_er <- rumour_base(
  graph = ger, n_iters = 100, inf_0 = 1,
  p_skep = 0.5, spr_rate = 0.5, rec_rate = 0.5,
  seed = FALSE, display = TRUE
)
res_cph <- rumour_base(
  graph = gcph, n_iters = 1000, inf_0 = 3,
  p_skep = 0.05, spr_rate = 0.3, rec_rate = 0.1, #thresh = 2,
  seed = FALSE, display = TRUE
)
gf <- make_full_graph(n = vcount(gcph))
res_f <- rumour_base(
  graph = gf, n_iters = 100, inf_0 = 1,
  p_skep = 0.2, spr_rate = 0.8, rec_rate = 0.1,
  seed = FALSE, display = TRUE
)

k <- degree(gcph, mode = "out") |> unname()
cls <- c("sus", "rec")
k_sr <- rbindlist(lapply(cls, function(x) {
    data.frame(class = x, k = k[res_cph[[x]]])
}))[, class := factor(class, cls)]

ggplot(k_sr, aes(k, fill = class)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1, alpha = 0.5,
    position = "identity", boundary = 0
  )

# g <- make_graph("Zachary")
g <- gcph

n_sim <- 1000
nodes <- V(g)
sims <- matrix(0, ncol = n_sim, nrow = vcount(g))
for (i in seq_len(n_sim)) {
  res <- rumour_base(
    graph = g, n_iters = 100, inf_0 = 1,
    p_skep = 0.2, spr_rate = 0.8, rec_rate = 0.2,
    seed = FALSE, display = FALSE
  )

  sims[nodes[res$sus], i] <- "s"
  sims[nodes[res$inf], i] <- "i"
  sims[nodes[res$rec], i] <- "r"
}

sir_freqs <- rbindlist(
  apply(
    sims, 1,
    \(r) list(sus = sum(r == "s"), inf = sum(r == "i"), rec = sum(r == "r"))
  )
)[, lapply(.SD, \(x) x / n_sim)]

# colour by frequency of final state = sus during n_sim runs
n_cols <- 20
freqs <- sir_freqs$sus
hist(freqs)
cols <- viridis::plasma(n_cols + 1)[
  1 + floor(n_cols * (freqs - min(freqs)) / (max(freqs) - min(freqs)))]
plot(
  g, layout = layout_with_graphopt(g),
  vertex.color = cols, vertex.size = 2 * degree(g)^0.3,
  vertex.label = NA, vertex.frame.color = NA
)

# colour by final state in sample run
res <- rumour_base(
  graph = g, n_iters = 100, inf_0 = 1,
  p_skep = 0.2, spr_rate = 0.8, rec_rate = 0.5,
  seed = FALSE, display = TRUE
)
V(g)$state <- ifelse(res$sus, "sus", ifelse(res$rec, "rec", "inf"))
cols <- rep("darkgoldenrod1", vcount(gcph))
cols[res$inf] <- "firebrick"
cols[res$rec] <- "dodgerblue4"
cols[res$start] <- "black"
plot(
  g, layout = layout_with_kk(g),
  vertex.color = cols, vertex.size = 2 * sqrt(degree(g)),
  vertex.label = NA, vertex.frame.color = NA
)
