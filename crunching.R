library(igraph)
library(data.table)
library(parallel)
source("src/utils.R")
source("src/rumour_base.R")
source("src/rumour_dose.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

calls <- fread(
  "data/copenhagen/calls.csv",
  col.names = c("time", "from", "to", "duration")
)[duration > 0, .N, keyby = .(a = pmin(from, to), b = pmax(from, to))]

sms <- fread(
  "data/copenhagen/sms.csv",
  col.names = c("time", "from", "to")
)[, .N, keyby = .(a = pmin(from, to), b = pmax(from, to))]

cph <- merge(calls, sms, all = TRUE)[
  is.na(N.x), N.x := 0][is.na(N.y), N.y := 0][
    , .(a, b, weight = 1 + log(N.x + N.y))][
      , weight := weight / max(weight)]

g <- graph_from_data_frame(cph, directed = FALSE) |>
  get_lcc() |>
  simplify()
V(g)$name <- seq_len(vcount(g))

get_results <- function(spr, rec, func, n_sim = 100) {
  res <- lapply(
    seq_len(n_sim),
    function(i) {
      if (func == "base") {
        rumour_base(
          graph = g, n_iters = 1e4, inf_0 = 1,
          p_skep = 0, spr_rate = spr, rec_rate = rec,
          seed = FALSE, display = FALSE
        )
      } else {
        rumour_dose(
          graph = g, n_iters = 1e4, inf_0 = 1,
          p_skep = 0, spr_rate = spr, rec_rate = rec, thresh = 2,
          seed = FALSE, display = FALSE
        )
      }
    }
  ) |>
    Reduce(
      function(x, y) {
        list(
          att_rate = x$att_rate + sum(y$reached) / n_sim,
          max_inf = x$max_inf + max(y$n_inf) / n_sim,
          duration = x$duration + y$duration / n_sim,
          reached = x$reached + 1 * y$reached / n_sim
        )
      },
      x = _,
      init = list(att_rate = 0, max_inf = 0, duration = 0, reached = 0)
    )

  res$spr_rate <- spr
  res$rec_rate <- rec

  return(res)
}

###################
# PARALLELIZATION #
###################

detectCores()
n_cores <- 4

spr_rates <- c(0.3, 0.7, 0.9)
rec_rates <- c(0.1, 0.2, 0.3)
rates <- expand.grid(spr = spr_rates, rec = rec_rates) |>
  transpose()
setnames(rates, as.character(1:ncol(rates)))

if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  results <- mclapply(
    rates, \(x) get_results(x[1], x[2], func = "dose", n_sim = 5),
    mc.cores = n_cores
  )
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster, c("get_results", "rumour_base", "rumour_dose", "rates", "g")
  )
  results <- parLapplyLB(
    cluster, rates,
    \(x) get_results(x[1], x[2], func = "dose", n_sim = 5)
  )
  stopCluster(cluster)
}

str(results)
dput(results, "simulations/dose.txt")

##################
# SERIAL VERSION #
##################

spr_rates <- c(0.3, 0.7, 0.9)
rec_rates <- c(0.1, 0.2, 0.3)

sttime <- Sys.time()
results <- outer(
  # spr_rate vector
  spr_rates,
  # rec_rate vector
  rec_rates,
  Vectorize(
    \(s, r) get_results(s, r, func = "dose", n_sim = 5),
    SIMPLIFY = FALSE
  )
)
dim(results) <- NULL # convert to list
fstime <- Sys.time()
(elapsed <- sttime - fstime)

# write to file
dput(results, "base.txt") # or "dose.txt"
# read back with
# res <- dget("out.txt")
