logbins <- function(data, base = 10, dlog = 0.1) {
  # determine counts from data
  values <- sort(unique(data))
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
