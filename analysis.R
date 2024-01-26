library(igraph)
library(ggplot2)
library(data.table)
source("src/utils.R")
source("src/rumour_base.R")
source("src/rumour_dose.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

g <- read_graph("data/graph_cph.graphml", format = "graphml")
k <- unname(degree(g))

base1 <- dget("out/base_011.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()
base5 <- dget("out/base_595.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()

base <- rbind(base1, base5)
cols <- names(base)[names(base) != "reached"]
base[, (cols) := lapply(.SD, unlist), .SDcols = cols]
base$att_rate <- base$att_rate / vcount(g)
setcolorder(base, c(5, 6, 1:4))
setorder(base, spr_rate, rec_rate)

fwrite(base, "out/base.csv")

dose1 <- dget("out/dose_011.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()
dose5 <- dget("out/dose_595.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()

dose <- rbind(dose1, dose5)
cols <- names(dose)[names(dose) != "reached"]
dose[, (cols) := lapply(.SD, unlist), .SDcols = cols]
dose$att_rate <- dose$att_rate / vcount(g)
setcolorder(dose, c(5, 6, 1:4))
setorder(dose, spr_rate, rec_rate)

fwrite(dose, "out/dose.csv")

############
# ANALYSIS #
############

base5 <- base[grep(".[0-9]5$", spr_rate)][grep(".[0-9]5$", rec_rate)] 
dose5 <- dose[grep(".[0-9]5$", spr_rate)][grep(".[0-9]5$", rec_rate)] 

base5 |>
  ggplot(aes(spr_rate, rec_rate, fill = att_rate)) +
    scale_fill_viridis_c() +
    geom_tile() +
    geom_line(aes(x = spr_rate, y = spr_rate))

dose5 |>
  ggplot(aes(spr_rate, rec_rate, fill = att_rate)) +
    scale_fill_viridis_c() +
    geom_tile()

base5 |>
  ggplot(aes(spr_rate, rec_rate, fill = log10(duration))) +
    scale_fill_viridis_c() +
    geom_tile()
    
base[rec_rate == 0.05] |>
  ggplot(aes(spr_rate, duration)) +
    geom_point()

dose5 |>
  ggplot(aes(spr_rate, rec_rate, fill = log10(duration))) +
    scale_fill_viridis_c() +
    geom_tile()

dose[rec_rate == 0.15] |>
  ggplot(aes(spr_rate, duration)) +
    geom_line()

base5 |>
  ggplot(aes(spr_rate, rec_rate, fill = max_inf)) +
    scale_fill_viridis_c() +
    geom_tile()

dose5 |>
  ggplot(aes(spr_rate, rec_rate, fill = max_inf)) +
    scale_fill_viridis_c() +
    geom_tile()


plot(k, dose[spr_rate == 0.25 & rec_rate == 0.05, 500 * unlist(reached)])

dose[spr_rate == 0.25 & rec_rate == 0.05,
     .(deg = degree(g), reached = unlist(reached))][
     , .(reached = median(reached), sd = sd(reached)), keyby = deg] |>
  ggplot(aes(deg, reached)) +
    geom_point()

dose[spr_rate == 0.25 & rec_rate == 0.05]

base5[att_rate > 0.1,
      .(cor = cor(unlist(reached), k)),
      keyby = .(spr_rate, rec_rate)] |>
  ggplot(aes(spr_rate, rec_rate, fill = cor)) +
    scale_fill_viridis_c() +
    geom_tile()

dose5[att_rate > 0.1, .(cor = cor(unlist(reached), k)), keyby = .(spr_rate, rec_rate)] |>
  ggplot(aes(spr_rate, rec_rate, fill = cor)) +
    scale_fill_viridis_c() +
    geom_tile()

base_cor <- base[att_rate > 0.1, .(cor = cor(unlist(reached), k), duration),
                 keyby = .(spr_rate, rec_rate)]
dose_cor <- dose[att_rate > 0.1, .(cor = cor(unlist(reached), k), duration),
                 keyby = .(spr_rate, rec_rate)]

base_bc <- base_cor[cor == max(cor), .(spr = spr_rate, rec = rec_rate)]
dose_bc <- dose_cor[cor == max(cor), .(spr = spr_rate, rec = rec_rate)]

res <- rumour_dose(
  graph = g, n_iters = 1e4, inf_0 = 1,
  p_skep = 0, spr_rate = 0.8, rec_rate = 0.2, thresh = 5,
  seed = FALSE, display = TRUE
)
V(g)$state <- ifelse(res$sus, "sus", "rec")
cols <- rep("darkgoldenrod1", vcount(g))
cols[res$rec] <- "dodgerblue4"
cols[res$start] <- "olivedrab3"
size <- 2 * sqrt(k)
size[res$start] <- 0.5 * max(size)
plot(
  g, layout = layout_as_tree(g, root = res$start),
  vertex.color = cols, vertex.size = 3,
  vertex.label = NA, vertex.frame.color = NA
)
