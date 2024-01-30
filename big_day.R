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
save_plot("deg_hist", img_dir, scale = 1.8)

# weights histogram
ggplot() +
  geom_histogram(
    aes(w, after_stat(density)),
    binwidth = 0.05, boundary = 0, fill = sir_pal[2],
    colour = bg_col
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = "Edge weight", y = "Density") +
  theme_sir()
save_plot("wgt_hist", img_dir, scale = 1.8)

# infection videos
base_vid <- rumour_base(
  graph = g, n_iters = 1e4, inf_0 = which.max(k),
  p_skep = 0.1, spr_rate = 0.85, rec_rate = 0.15,
  seed = 80085, display = TRUE, save_plots = TRUE,
  graph_lay = g_lay
)
system("ffmpeg -r 10 -i img/vid/base_%04d.jpg img/vid/base.mp4")
dose_vid <- rumour_dose(
  graph = g, n_iters = 1e4, inf_0 = which.max(k),
  p_skep = 0.1, spr_rate = 0.85, rec_rate = 0.15, thresh = 5,
  seed = 80085, display = TRUE, save_plots = TRUE,
  graph_lay = g_lay
)
system("ffmpeg -r 10 -i img/vid/dose_%04d.jpg img/vid/dose.mp4")

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

# attack rate heatmaps
min_ar <- min(base_grid[, min(att_rate)], dose_grid[, min(att_rate)])
max_ar <- max(base_grid[, max(att_rate)], dose_grid[, max(att_rate)])
ggplot(base_grid, aes(spr_rate, rec_rate, fill = att_rate)) +
  scale_fill_viridis_c(lim = c(min_ar, max_ar)) +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Final attack rate",
    title = "Simple contagion"
  ) +
  theme_sir()
save_plot("att_rate_grid_base", img_dir)
ggplot(dose_grid, aes(spr_rate, rec_rate, fill = att_rate)) +
  scale_fill_viridis_c(lim = c(min_ar, max_ar)) +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Final attack rate",
    title = "Complex contagion"
  ) +
  theme_sir()
save_plot("att_rate_grid_dose", img_dir)

# duration heatmaps
min_dur <- min(base_grid[, min(duration)], dose_grid[, min(duration)])
max_dur <- max(base_grid[, max(duration)], dose_grid[, max(duration)])
ggplot(base_grid, aes(spr_rate, rec_rate, fill = duration)) +
  scale_fill_viridis_c(lim = c(min_dur, max_dur)) +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Epidemic duration",
    title = "Simple contagion"
  ) +
  theme_sir()
save_plot("duration_grid_base", img_dir)
ggplot(dose_grid, aes(spr_rate, rec_rate, fill = duration)) +
  scale_fill_viridis_c() +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Epidemic duration",
    title = "Complex contagion"
  ) +
  theme_sir()
save_plot("duration_grid_dose", img_dir)

# max_inf heatmaps
min_mi <- min(base_grid[, min(max_inf)], dose_grid[, min(max_inf)])
max_mi <- max(base_grid[, max(max_inf)], dose_grid[, max(max_inf)])
ggplot(base_grid, aes(spr_rate, rec_rate, fill = max_inf)) +
  scale_fill_viridis_c(lim = c(min_mi, max_mi)) +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Peak prevalence",
    title = "Simple contagion"
  ) +
  theme_sir()
save_plot("max_inf_grid_base", img_dir)
ggplot(dose_grid, aes(spr_rate, rec_rate, fill = max_inf)) +
  scale_fill_viridis_c(lim = c(min_mi, max_mi)) +
  geom_tile() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Spreading rate",
    y = "Recovery rate",
    fill = "Peak prevalence",
    title = "Complex contagion"
  ) +
  theme_sir()
save_plot("max_inf_grid_dose", img_dir)

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

base_sk[, .(duration_mu = mean(duration), duration_sd = sd(duration),
            att_rate_mu = mean(att_rate), att_rate_sd = sd(att_rate),
            max_inf_mu = mean(max_inf), max_inf_sd = sd(max_inf)),
        by = psk] |>
  melt(id = "psk", measure = patterns(mean = "_mu$", sd = "_sd$")) |>
  ggplot(aes(psk, mean)) +
    geom_line() +
    labs(
      x = "Skeptical probability", y = NULL,
      title = "Simple contagion"
    ) +
    facet_grid(
      rows = vars(variable), scales = "free_y",
      labeller = as_labeller(
        c(`1` = "Epid. duration", `2` = "Attack rate", `3` = "Peak prevalence")
      )
    ) +
    theme_sir()

################
# NODE BY NODE #
################

nbn_files <- list.files("out/nbn", full.names = TRUE)
base_nbn <- purrr::map(
  nbn_files[grep("base", nbn_files)], dget,
  .progress = TRUE
) |>
  unlist(recursive = FALSE, use.names = FALSE)
dose_nbn <- purrr::map(
  nbn_files[grep("dose", nbn_files)], dget,
  .progress = TRUE
) |>
  unlist(recursive = FALSE, use.names = FALSE)

# final attack rate by starting node
base_ar <- rbindlist(
  lapply(base_nbn, \(x) list(start = x$start, att_rate = sum(x$reached)))
)[, .(att_rate = mean(att_rate) / n_nodes), keyby = start]
dose_ar <- rbindlist(
  lapply(dose_nbn, \(x) list(start = x$start, att_rate = sum(x$reached)))
)[, .(att_rate = mean(att_rate) / n_nodes), keyby = start]

# correlation between closeness and final attack rate
prank <- unname(page_rank(g)$vector)
close <- unname(closeness(g))
eigen <- unname(eigen_centrality(g)$vector)
data.table(att_rate = base_ar$att_rate, prank, close, eigen = log10(eigen)) |>
  melt(
    measure.vars = c("prank", "close", "eigen"),
    variable.name = "ctype", value.name = "cent"
  ) |>
  ggplot(aes(att_rate, cent)) +
    geom_point() +
    facet_grid(rows = vars(ctype), scales = "free_y")
data.table(att_rate = dose_ar$att_rate, prank, close, eigen = log10(eigen)) |>
  melt(
    measure.vars = c("prank", "close", "eigen"),
    variable.name = "ctype", value.name = "cent"
  ) |>
  ggplot(aes(att_rate, cent)) +
    geom_point() +
    facet_grid(rows = vars(ctype), scales = "free_y")

# node distance vs when_inf
dists <- distances(g)
base_wi <- rbindlist(
  lapply(
    base_nbn,
    \(x) list(start = x$start, end = 1:n_nodes, when = x$when_inf)
  )
)[when > 0, .(when = mean(when), .N), keyby = .(start, end)]
  , dist := mapply(\(a, b) dists[a, b], start, end)]
ggplot(base_wi, aes(dist, N)) +
  geom_point(alpha = 0.1, size = sz_micro)

base_wi[, .(when = mean(when)), keyby = .(k1 = k[start], k2 = k[end])] |>
  ggplot(aes(k1, k2, fill = when)) +
    scale_fill_viridis_c() +
    geom_tile() +
    coord_cartesian(expand = FALSE) +
    theme_sir()

ggplot(base_wi, aes(start, end, fill = when)) +
  scale_fill_viridis_c()
  geom_tile() +

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

#############
# EVOLUTION #
#############

base <- rumour_base(
  graph = g, n_iters = 1e4, inf_0 = 0,
  p_skep = 0.1, spr_rate = 0.85, rec_rate = 0.15,
  display = TRUE
)
