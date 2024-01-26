library(igraph)
library(data.table)

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
setcolorder(dose, c(5, 6, 1:4))
setorder(dose, spr_rate, rec_rate)

fwrite(dose, "out/dose.csv")
