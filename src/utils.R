library(showtext)
font_add_google("Jost", "plot-font")
showtext_auto()

# palette
sir_pal <- c("#f2e9e4", "#7180b9", "#22223b", "#a50104", "#058a5e")
bg_col <- sir_pal[1]
fg_col <- sir_pal[3]

# for plotting labels
nice_sir <- c("Susceptible", "Infected", "Recovered")

# text/line sizes
base_pt <- 15
sz_big <- ggplot2::rel(1.778)
sz_medium <- ggplot2::rel(1.333)
sz_normal <- ggplot2::rel(1)
sz_small <- ggplot2::rel(0.75)
sz_tiny <- ggplot2::rel(0.563)
sz_micro <- ggplot2::rel(0.422)

theme_sir <- function() {
  library(ggplot2)

  theme_classic(
    base_family = "plot-font",
    base_size = base_pt,
  ) +
  theme(
    plot.background = element_rect(fill = bg_col, colour = NA),
    panel.background = element_rect(fill = bg_col, colour = NA),
    legend.background = element_rect(fill = bg_col, colour = NA),
    # main title
    plot.title = element_text(
      colour = fg_col,
      size = sz_big,
      margin = margin(10, 0, 10, 0),
    ),
    plot.subtitle = element_text(
      colour = fg_col,
      size = sz_normal,
      margin = margin(0, 0, 20, 0),
    ),
    # axes titles and labels
    axis.title = element_text(
      size = sz_normal,
      hjust = 1,
      colour = fg_col,
    ),
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text = element_text(size = sz_small, colour = fg_col),
    # axis lines and ticks
    axis.line = element_line(linewidth = sz_tiny, colour = fg_col),
    axis.ticks = element_line(linewidth = sz_tiny, colour = fg_col),
    # legend titles and labels
    legend.title = element_text(
      size = sz_normal,
      colour = fg_col,
    ),
    legend.text = element_text(size = sz_small, colour = fg_col),
    # caption
    plot.caption = element_text(size = sz_small, hjust = 0, colour = fg_col),
    # facet labels
    strip.background = element_rect(fill = bg_col, colour = NA),
    strip.text = element_text(size = sz_normal, colour = fg_col)
  )
}

# common image saving
save_plot <- function(
  fname = "plot",
  prefix = getwd(),
  plot = ggplot2::last_plot(),
  asp_ratio = c(4, 3),
  scale = 2.5
) {
  ggplot2::ggsave(
    filename = sprintf("%s/%s.svg", prefix, fname),
    plot = plot,
    width = asp_ratio[1] * scale,
    height = asp_ratio[2] * scale,
    units = "in"
  )
}

ggsir <- function(n_inf, n_rec, n_sus) {
  data.frame(inf = n_inf, rec = n_rec, sus = n_sus) |>
    dplyr::mutate(iter = dplyr::row_number()) |>
    tidyr::pivot_longer(-"iter", names_to = "class", values_to = "pop") |>
    dplyr::mutate(class = factor(class, levels = c("sus", "inf", "rec"))) |>
    ggplot2::ggplot(ggplot2::aes(iter, pop, colour = class)) +
      ggplot2::geom_line(linewidth = sz_small) +
      ggplot2::scale_colour_manual(
        values = sir_pal[3:5],
        labels = c("Susceptible", "Infected", "Recovered")
      ) +
      ggplot2::labs(
        x = "Time step",
        y = "Fractional prevalence",
        colour = "Compartment"
      ) +
      theme_sir()
}

get_lcc <- function(graph, mode = "weak") {
  # extract the largest connected component from an igraph network
  components <- igraph::clusters(graph, mode = mode)
  lcc_id <- which.max(components$csize)
  lcc_verts <- igraph::V(graph)[components$membership == lcc_id]

  return(igraph::induced_subgraph(graph, lcc_verts))
}

logbins <- function(data, base = 10, dlog = 0.1) {
  # determine counts from data
  values <- seq_len(max(data))
  counts <- vapply(values, \(x) sum(data == x), integer(1))

  # get (log-)midpoints and lower/upper limits (to nearest integer)
  log_mid <- seq(0, log(max(values), base) + dlog, by = dlog)
  lower <- ceiling(base^(log_mid - 0.5 * dlog))
  upper <- floor(base^(log_mid + 0.5 * dlog))

  # select only bins larger than two, e.g. lower = 9, upper = 11
  bin_window <- which(upper > lower + 1)
  lower <- lower[bin_window]
  upper <- upper[bin_window]

  # take values before bin_window directly from counts
  x <- c(log(seq(1, lower[1] - 1), base), log_mid[bin_window])
  y <- c(
    counts[values < lower[1]],  # these are from counts
    # here we average the counts falling inside each remaining bin
    mapply(
      \(l, u) sum(counts[l <= values & values <= u]) / (1 + u - l),
      lower, upper
    )
  )

  y <- log(y, base)

  # remove rows with < 0 log-counts when returning
  return(data.frame(x = x[y > 0], y = y[y > 0]))
}
