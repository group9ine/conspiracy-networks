library(igraph)
library(data.table)
library(ggplot2)
source("src/utils.R")
source("src/rumour_dose.R")
source("src/rumour_base.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

img_dir <- "img"

#####################
# THIS IS A NETWORK #
#####################

g <- read_graph("data/graph_cph.graphml", format = "graphml")
k <- unname(degree(g))
w <- E(g)$weight
n_nodes <- vcount(g)
g_lay <- layout_with_graphopt(
  g, niter = 800, charge = 0.008, mass = 40, spring.length = 2
)

# colour nodes by centrality
cent <- eigen_centrality(g)$vector
n_cols <- 10
vcols <- colorRampPalette(sir_pal[c(3, 4)])(n_cols)[cut(
  cent, breaks = seq(min(cent), max(cent), len = n_cols + 1), include.lowest = TRUE
)]

png(paste(img_dir, "nw_plot.png", sep = "/"), 1000, 1000, bg = sir_pal[1])
plot(
  g, layout = g_lay,
  vertex.color = vcols, vertex.size = 2.25 * k^0.3,
  vertex.frame.color = NA, vertex.label = NA,
  edge.color = sir_pal[2], edge.width = 0.5 * w / min(w)
)
dev.off()

# degree histogram
ggplot() +
  geom_histogram(
    aes(k, after_stat(density)),
    binwidth = 1, boundary = 0.5, fill = sir_pal[2],
    colour = bg_col
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = "Degree", y = "Density") +
  theme_sir()
save_plot("deg_hist", img_dir)

# weights histogram
ggplot() +
  geom_histogram(
    aes(w, after_stat(density)),
    binwidth = 0.05, boundary = 0, fill = sir_pal[2],
    colour = bg_col
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),


#################
# GRID SEARCHES #
#################

# read block of simulation results from varying spr and rec
base_grid <- dget("out/base_595.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()

# unlist all columns but reached
unl <- names(base_grid)[names(base_grid) != "reached"]
base_grid[, (unl) := lapply(.SD, unlist), .SDcols = unl]
# normalize the att_rate
base_grid$att_rate <- base_grid$att_rate / n_nodes
# just some reordering
setcolorder(base_grid, c(5, 6, 1:4))
setorder(base_grid, spr_rate, rec_rate)

fwrite(base_grid, "out/base_grid.csv")

# do the same for dose
dose_grid <- dget("out/dose_595.txt") |>
  lapply(t) |>
  do.call(rbind, args = _) |>
  as.data.table()
unl <- names(dose_grid)[names(dose_grid) != "reached"]
dose_grid[, (unl) := lapply(.SD, unlist), .SDcols = unl]
dose_grid$att_rate <- dose_grid$att_rate / n_nodes
setcolorder(dose_grid, c(5, 6, 1:4))
setorder(dose_grid, spr_rate, rec_rate)
fwrite(dose_grid, "out/dose_grid.csv")

##### TODO analysis

############
# SKEPTICS #
############

base_sk <- fread("out/base_skeptics.csv")
str(base_sk)  # list columns are weird, need to convert them from char
atoi <- function(x) {
  as.integer(unlist(strsplit(x, "|", fixed = TRUE), FALSE, FALSE))
}

conv <- c("reached", "dir_rec", "when_inf")
base_sk[, (conv) := lapply(.SD, \(x) lapply(x, atoi)), .SDcols = conv][
  , `:=`(att_rate = att_rate / n_nodes, n_dir = sapply(dir_rec, sum) / n_nodes)]

# do the same for dose
dose_sk <- fread("out/dose_skeptics.csv")[
  , (conv) := lapply(.SD, \(x) lapply(x, atoi)), .SDcols = conv][
  , `:=`(att_rate = att_rate / n_nodes, n_dir = sapply(dir_rec, sum) / n_nodes)]

#### TODO analysis

#####################
# HOMOGENOUS MIXING #
#####################

hm_data <- fread("hma/arrays_data.csv")[, iter := .I]
means <- melt(
  hm_data, id.vars = "iter", measure.vars = c("s", "i", "r"),
  variable.name = "class", value.name = "mean"
)
sigmas <- melt(
  hm_data, id.vars = "iter", measure.vars = patterns("sigma"),
  variable.name = "class", value.name = "sigma"
)[, class := sub("sigma_", "", class, fixed = TRUE)]

hm_data <- merge(means, sigmas)[
  , class := factor(class, c("s", "i", "r"))][
  order(iter, class)]

ggplot(hm_data, aes(iter, mean,fill = class)) +
  geom_ribbon(aes(ymin = mean - sigma, ymax = mean + sigma), alpha = 0.3) +
  geom_line(aes(colour = class), linewidth = sz_small) +
  scale_colour_manual(values = sir_pal[3:5], labels = nice_sir) +
  scale_fill_manual(values = sir_pal[3:5], labels = nice_sir) +
  labs(
    x = "Iteration", y = "Fractional prevalence",
    colour = "Compartment", fill = "Compartment"
  ) +
  theme_sir()

save_plot("homo_mix", img_dir)

###############
# OTHER STUFF #
###############

rumour_base(
  graph = g, n_iters = 1e4, inf_0 = 0,
  p_skep = 0.1, spr_rate = 0.85, rec_rate = 0.15,
  display = TRUE
)
