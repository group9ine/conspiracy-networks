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

ggsir <- function(n_inf, n_rec, n_sus) {
  data.frame(inf = n_inf, rec = n_rec, sus = n_sus) |>
    dplyr::mutate(iter = dplyr::row_number()) |>
    tidyr::pivot_longer(-"iter", names_to = "class", values_to = "pop") |>
    dplyr::mutate(class = factor(class, levels = c("sus", "inf", "rec"))) |>
    ggplot2::ggplot(ggplot2::aes(iter, pop, colour = class)) +
      ggplot2::geom_line() +
      ggplot2::scale_colour_manual(
        values = c("darkgoldenrod1", "firebrick", "dodgerblue4"),
        labels = c("Susceptible", "Infected", "Recovered")
      ) +
      ggplot2::labs(
        x = "Time step",
        y = "Fractional prevalence",
        colour = "Compartment"
      )
}

get_lcc <- function(graph, mode = "weak") {
  # extract the largest connected component from an igraph network
  components <- igraph::clusters(graph, mode = mode)
  lcc_id <- which.max(components$csize)
  lcc_verts <- igraph::V(graph)[components$membership == lcc_id]

  return(igraph::induced_subgraph(graph, lcc_verts))
}
