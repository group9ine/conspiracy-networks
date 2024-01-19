library(igraph)
library(ggplot2)
library(data.table)
source("src/utils.R")
source("src/rumours_rnd.R")
source("src/rumours_skep.R")

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
dcph <- merge(dc, merge(ds, df, all = TRUE), all = TRUE)[
  is.na(N.x), N.x := 0][is.na(N.y), N.y := 0][
      N.x + N.y > 0, weight := 1 + log(N.x + N.y)][
        N.x + N.y == 0, weight := 1][
          , c("N.x", "N.y") := NULL]

gcph <- graph_from_data_frame(dcph, directed = FALSE) |> simplify()
V(gcph)$name <- seq_len(vcount(gcph))  # normalize naming
is_simple(gcph)
clusters(gcph)  # only one cc, good

# null model
ger <- sample_gnm(n = vcount(gcph), m = ecount(gcph))
res_er <- rumour_skep(
  graph = ger, n_iters = 100, inf_0 = 1,
  p_talk = 0.8, p_skep = 0.2, p_stop = 0.2,
  d_wind = 7, thresh = 2,
  seed = FALSE, display = TRUE
)
res_cph <- rumour_skep(
  graph = gcph, n_iters = 100, inf_0 = 1,
  p_talk = 0.6, p_skep = 0.2, p_stop = 0.2,
  d_wind = 7, thresh = 2,
  seed = FALSE, display = TRUE
)

k <- degree(gcph, mode = "out") |> unname()
k_sr <- rbindlist(
  lapply(c("sus", "rec"), \(x) data.frame(class = x, k = k[res_cph[[x]]]))
)[, class := factor(class, c("sus", "rec"))]

ggplot(k_sr, aes(k, fill = class)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1, alpha = 0.5,
    position = "identity", boundary = 0
  )
