library(igraph)
library(data.table)

base1 <- dget("out/base_011.txt")
base5 <- dget("out/base_595.txt")

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
setcolorder(base, c(5, 6, 1:4))
setorder(base, spr_rate, rec_rate)
