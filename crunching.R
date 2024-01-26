library(igraph)
library(data.table)
library(parallel)
source("src/utils.R")
source("src/rumour_base.R")
source("src/rumour_dose.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

g <- read_graph("data/graph_cph.graphml", format = "graphml")

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
          p_skep = 0, spr_rate = spr, rec_rate = rec, thresh = 5,
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
n_cores <- 3

spr_rates <- seq(0.1, 1, 0.1)
rec_rates <- seq(0.1, 1, 0.1)
rates <- expand.grid(spr = spr_rates, rec = rec_rates) |>
  transpose()
setnames(rates, as.character(1:ncol(rates)))

if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  st_time <- Sys.time()
  results <- mclapply(
    rates, \(x) get_results(x[1], x[2], func = "dose", n_sim = 500),
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  fs_time <- Sys.time()
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster, c("get_results", "rumour_base", "rumour_dose", "rates", "g")
  )
  st_time <- Sys.time()
  results <- parLapplyLB(
    cluster, rates,
    \(x) get_results(x[1], x[2], func = "dose", n_sim = 5)
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}

fs_time - st_time

str(results)
dput(results, "simulations/dose.txt")


##################
# SERIAL VERSION #
##################

spr_rates <- c(0.5, 0.6)
rec_rates <- c(0.5, 0.6)

st_time <- Sys.time()
results <- outer(
  # spr_rate vector
  spr_rates,
  # rec_rate vector
  rec_rates,
  Vectorize(
    \(s, r) get_results(s, r, func = "dose", n_sim = 500),
    SIMPLIFY = FALSE
  )
)
dim(results) <- NULL # convert to list
fs_time <- Sys.time()
fs_time - st_time

# write to file
dput(results, "base.txt") # or "dose.txt"
# read back with
# res <- dget("out.txt")

#######################
# SKEPTICISM ANALYSIS #
#######################
