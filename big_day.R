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
prank <- page_rank(g)$vector
close <- closeness(g)
eigen <- eigen_centrality(g)$vector

n_cols <- 10
png(paste(img_dir, "nw_prank.png", sep = "/"), 1000, 1000, bg = sir_pal[1])
vcols <- colorRampPalette(sir_pal[c(3, 4)])(n_cols)[cut(
  prank, breaks = seq(min(prank), max(prank), len = n_cols + 1),
  include.lowest = TRUE
)]
plot(
  g, layout = g_lay,
  vertex.color = vcols, vertex.size = 2.25 * k^0.3,
  vertex.frame.color = NA, vertex.label = NA,
  edge.color = sir_pal[2], edge.width = 0.5 * w / min(w)
)
dev.off()
png(paste(img_dir, "nw_close.png", sep = "/"), 1000, 1000, bg = sir_pal[1])
vcols <- colorRampPalette(sir_pal[c(3, 4)])(n_cols)[cut(
  close, breaks = seq(min(close), max(close), len = n_cols + 1),
  include.lowest = TRUE
)]
plot(
  g, layout = g_lay,
  vertex.color = vcols, vertex.size = 2.25 * k^0.3,
  vertex.frame.color = NA, vertex.label = NA,
  edge.color = sir_pal[2], edge.width = 0.5 * w / min(w)
)
dev.off()
png(paste(img_dir, "nw_eigen.png", sep = "/"), 1000, 1000, bg = sir_pal[1])
vcols <- colorRampPalette(sir_pal[c(3, 4)])(n_cols)[cut(
  eigen, breaks = seq(min(eigen), max(eigen), len = n_cols + 1),
  include.lowest = TRUE
)]
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

base_sk[, .(variable = names(.SD),
            mean = sapply(.SD, mean),
            sd = sapply(.SD, sd)),
        by = psk, .SDcols = c("duration", "att_rate", "max_inf")] |>
  ggplot(aes(psk, mean)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      fill = fg_col, alpha = 0.3
    ) +
    geom_line(colour = fg_col, linewidth = sz_small) +
    labs(
      x = "Skeptical probability", y = NULL,
      title = "Simple contagion"
    ) +
    facet_grid(
      rows = vars(variable), scales = "free_y",
      labeller = as_labeller(
        c(duration = "Epidemic duration",
          att_rate = "Final attack rate",
          max_inf = "Peak prevalence")
      )
    ) +
    theme_sir()
save_plot("skep_base", img_dir, asp_ratio = c(3, 4))
dose_sk[, .(variable = names(.SD),
            mean = sapply(.SD, mean),
            sd = sapply(.SD, sd)),
        by = psk, .SDcols = c("duration", "att_rate", "max_inf")] |>
  ggplot(aes(psk, mean)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      fill = fg_col, alpha = 0.3
    ) +
    geom_line(colour = fg_col, linewidth = sz_small) +
    labs(
      x = "Skeptical probability", y = NULL,
      title = "Complex contagion"
    ) +
    facet_grid(
      rows = vars(variable), scales = "free_y",
      labeller = as_labeller(
        c(duration = "Epidemic duration",
          att_rate = "Final attack rate",
          max_inf = "Peak prevalence")
      )
    ) +
    theme_sir()
save_plot("skep_dose", img_dir, asp_ratio = c(3, 4))

#################
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

base_arcen <- data.table(
  att_rate = base_ar$att_rate, prank, close, eigen = log10(eigen)
)
# print correlations
base_arcen[, lapply(.SD, \(x) cor(att_rate, x)), .SDcols = !"att_rate"]
# plot
melt(
  base_arcen,
  measure.vars = c("prank", "close", "eigen"),
  variable.name = "ctype", value.name = "cent"
) |>
  ggplot(aes(att_rate, cent)) +
    geom_point(colour = fg_col, size = sz_small) +
    labs(
      x = "Final attack rate", y = "Centrality",
      title = "Simple contagion"
    ) +
    facet_grid(
      rows = vars(ctype), scales = "free_y",
      labeller = as_labeller(c(
        prank = "PageRank", close = "Closeness",
        eigen = "Eigenvector (log)"
      ))
    ) +
    theme_sir()
save_plot("cent_arate_base", img_dir, asp_ratio = c(3, 4))

dose_arcen <- data.table(
  att_rate = dose_ar$att_rate, prank, close, eigen = log10(eigen)
)
# print correlations
dose_arcen[, lapply(.SD, \(x) cor(att_rate, x)), .SDcols = !"att_rate"]
# plot
melt(
  dose_arcen,
  measure.vars = c("prank", "close", "eigen"),
  variable.name = "ctype", value.name = "cent"
) |>
  ggplot(aes(att_rate, cent)) +
    geom_point(colour = fg_col, size = sz_small) +
    labs(
      x = "Final attack rate", y = "Centrality",
      title = "Complex contagion"
    ) +
    facet_grid(
      rows = vars(ctype), scales = "free_y",
      labeller = as_labeller(c(
        prank = "PageRank", close = "Closeness",
        eigen = "Eigenvector (log)"
      ))
    ) +
    theme_sir()
save_plot("cent_arate_dose", img_dir, asp_ratio = c(3, 4))

# node distance vs when_inf
dists <- distances(g)
base_wi <- rbindlist(
  lapply(
    base_nbn,
    \(x) list(start = x$start, end = 1:n_nodes, when = x$when_inf)
  ) 
)[when > 0, .(when = mean(when), .N), keyby = .(start, end)][
  , dist := mapply(\(a, b) dists[a, b], start, end)]
ggplot(base_wi, aes(dist, when)) +
  geom_bin2d(aes(fill = after_stat(density)), bins = 50) +
  scale_fill_viridis_c() +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Weighted node distance",
    y = "Average time to infection",
    title = "Simple contagion",
    fill = "Density"
  ) +
  theme_sir()
save_plot("tinf_dist_base", img_dir)

dose_wi <- rbindlist(
  lapply(
    dose_nbn,
    \(x) list(start = x$start, end = 1:n_nodes, when = x$when_inf)
  ) 
)[when > 0, .(when = mean(when), .N), keyby = .(start, end)][
  , dist := mapply(\(a, b) dists[a, b], start, end)]
ggplot(dose_wi, aes(dist, when)) +
  geom_bin2d(aes(fill = after_stat(density)), bins = 50) +
  scale_fill_viridis_c() +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Weighted node distance",
    y = "Average time to infection",
    title = "Complex contagion",
    fill = "Density"
  ) +
  theme_sir()
save_plot("tinf_dist_dose", img_dir)

base_wi[, .(when = mean(when)), keyby = .(k1 = k[start], k2 = k[end])] |>
  ggplot(aes(k1, k2, fill = when)) +
    scale_fill_viridis_c(breaks = seq(5, max(base_wi$when), 10)) +
    geom_tile() +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "Degree of starting node",
      y = "Degree of infected node",
      fill = "Avg. time to inf.",
      title = "Simple contagion"
    ) +
    theme_sir()
save_plot("tinf_degs_base", img_dir)
dose_wi[, .(when = mean(when)), keyby = .(k1 = k[start], k2 = k[end])] |>
  ggplot(aes(k1, k2, fill = when)) +
    scale_fill_viridis_c(breaks = seq(50, max(dose_wi$when), 50)) +
    geom_tile() +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "Degree of starting node",
      y = "Degree of infected node",
      fill = "Avg. time to inf.",
      title = "Simple contagion"
    ) +
    theme_sir()
save_plot("tinf_degs_dose", img_dir)

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
    x = "Iteration", y = "Fractional occupation",
    colour = "Compartment", fill = "Compartment",
    title = "Stochastic homogeneous mixing"
  ) +
  theme_sir()
save_plot("homo_mix", img_dir)

#############
# EVOLUTION #
#############

base_evo <- lapply(seq_along(base_nbn), function(i) {
  x <- base_nbn[[i]]
  if (sum(x$reached) > 0.2 * n_nodes) {
    return(list(
      sim = i, iter = seq_along(x$n_inf),
      sus = x$n_sus, inf = x$n_inf, rec = x$n_rec
    ))
  }
}) |>
  rbindlist() |>
  melt(id = c("iter", "sim"), variable.name = "class")

base_evo[sim %in% sample(sim, 500) & iter < 200][
  , med := median(value), by = .(iter, class)] |>
  ggplot(aes(iter, value, colour = class)) +
    geom_line(
      aes(group = interaction(sim, class)),
      linewidth = sz_micro, alpha = 0.05
    ) +
    geom_line(aes(iter, med), linewidth = sz_normal) +
    scale_colour_manual(values = sir_pal[3:5], labels = nice_sir) +
    labs(
      x = "Iteration", y = "Fractional occupation",
      colour = "Compartment", title = "Simple contagion",
      subtitle = "500 stochastic trajectories"
    ) +
    theme_sir()
save_plot("ghost_base", img_dir)

dose_evo <- lapply(seq_along(dose_nbn), function(i) {
  x <- dose_nbn[[i]]
  if (sum(x$reached) > 0.05 * n_nodes) {
    return(list(
      sim = i, iter = seq_along(x$n_inf),
      sus = x$n_sus, inf = x$n_inf, rec = x$n_rec
    ))
  }
}) |>
  rbindlist() |>
  melt(id = c("iter", "sim"), variable.name = "class")

dose_evo[sim %in% sample(sim, 500) & iter < 800][
  , med := median(value), by = .(iter, class)] |>
  ggplot(aes(iter, value, colour = class)) +
    geom_line(
      aes(group = interaction(sim, class)),
      linewidth = sz_micro, alpha = 0.05
    ) +
    geom_line(aes(iter, med), linewidth = sz_normal) +
    scale_colour_manual(values = sir_pal[3:5], labels = nice_sir) +
    labs(
      x = "Iteration", y = "Fractional occupation",
      colour = "Compartment", title = "Complex contagion",
      subtitle = "500 stochastic trajectories"
    ) +
    theme_sir()
save_plot("ghost_dose", img_dir)

base_kinf <- lapply(seq_along(base_nbn), function(i) {
  x <- base_nbn[[i]]
  if (sum(x$reached) > 0.2 * n_nodes) {
    k_inf <- x$k_inf[!is.na(x$k_inf)]
    return(list(sim = i, iter = seq_along(k_inf), k = k_inf))
  }
}) |>
  rbindlist()

base_kinf[sim %in% sample(sim, 750) & iter < 150][
  , med := median(k), by = iter] |>
  ggplot(aes(iter, k)) +
    geom_line(
      aes(group = sim),
      linewidth = sz_micro, alpha = 0.07,
      colour = sir_pal[2]
    ) +
    geom_line(
      aes(y = med), linewidth = sz_normal,
      colour = sir_pal[3]
    ) +
    labs(
      x = "Iteration", y = "Average degree of newly infected",
      title = "Simple contagion", subtitle = "750 stochastic trajectories"
    ) +
    theme_sir()
save_plot("kinf_base", img_dir)

dose_kinf <- lapply(seq_along(dose_nbn), function(i) {
  x <- dose_nbn[[i]]
  if (sum(x$reached) > 0.05 * n_nodes) {
    k_inf <- x$k_inf[!is.na(x$k_inf)]
    return(list(sim = i, iter = seq_along(k_inf), k = k_inf))
  }
}) |>
  rbindlist()

dose_kinf[sim %in% sample(sim, 750) & iter < 500][
  , med := median(k), by = iter] |>
  ggplot(aes(iter, k)) +
    geom_line(
      aes(group = sim),
      linewidth = sz_micro, alpha = 0.07,
      colour = sir_pal[2]
    ) +
    geom_line(
      aes(y = med), linewidth = sz_normal,
      colour = sir_pal[3]
    ) +
    labs(
      x = "Iteration", y = "Average degree of newly infected",
      title = "Complex contagion", subtitle = "750 stochastic trajectories"
    ) +
    theme_sir()
save_plot("kinf_dose", img_dir)

#########################
# NULL MODEL COMPARISON #
#########################

nbn_files <- list.files("out/nbn", full.names = TRUE)
null_nbn <- purrr::map(
  nbn_files[grep("null", nbn_files)], dget,
  .progress = TRUE
) |>
  unlist(recursive = FALSE, use.names = FALSE)

null_ar <- lapply(null_nbn, \(x) sum(x$reached) / n_nodes) |>
  unlist(FALSE, FALSE)

base_ar <- lapply(base_nbn, \(x) sum(x$reached) / n_nodes) |>
  unlist(FALSE, FALSE)

data.table(null = null_ar, base = base_ar[seq_along(null_ar)]) |>
  melt(measure = c("null", "base")) |>
  ggplot(aes(x = value, fill = variable)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = 0.025, boundary = 0, alpha = 0.3,
      position = "identity"
    ) +
    scale_fill_manual(
      values = c(null = sir_pal[3], base = sir_pal[4]),
      labels = c(null = "Null model", base = "Copenhagen")
    ) +
    labs(x = "Final attack rate", y = "Density", fill = NULL, title = "Final attack rate histograms", subtitle = "Copenhagen network vs its configuration null model") +
    coord_cartesian(expand = FALSE) +
    theme_sir()
save_plot("null_model", "img")

ggplot() +
  geom_histogram(
    aes(x = null_ar, y = after_stat(density)),
    binwidth = 0.025, boundary = 0
  )

##########
# EXTRAS #
##########

rinf <- function(x, a, b, p){
  return(a/b*(1-x) - a/b*log(1-x) + p - a/b - p*(1-x) - x)
}

x <- seq(0, by = 0.001, len = 1000)

data.table(x = x, p0 = rinf(x, 0.15, 0.3, 0.1), p1 = rinf(x, 0.15, 0.3, 1)) |>
  melt(id = "x") |>
  ggplot(aes(x, value, colour = variable)) +
    geom_hline(
      aes(yintercept = 0), linewidth = sz_small,
      linetype = "dashed", colour = sir_pal[2]
    ) +
    geom_line(linewidth = sz_small) +
    scale_colour_manual(
      values = c(p0 = sir_pal[3], p1 = sir_pal[4]),
      labels = c(p0 = expression(paste(italic("p"), " = 0.1")),
                 p1 = expression(paste(italic("p"), " = 1.0")))
    ) +
    labs(
      x = expression(italic("r")),
      y = expression(paste(italic("f"), "(", italic("r"), ")")),
      colour = NULL
    ) +
    theme_sir()
save_plot("universality", img_dir, asp_ratio = c(4.5, 3))
