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

get_results <- function(start, psk, spr, rec, func, n_sim = 100) {
  res <- lapply(
    seq_len(n_sim),
    function(i) {
      if (func == "base") {
        rumour_base(
          graph = g, n_iters = 1e4, inf_0 = start,
          p_skep = psk, spr_rate = spr, rec_rate = rec,
          seed = FALSE, display = FALSE
        )
      } else {
        rumour_dose(
          graph = g, n_iters = 1e4, inf_0 = start,
          p_skep = psk, spr_rate = spr, rec_rate = rec, thresh = 5,
          seed = FALSE, display = FALSE
        )
      }
    }
  )
  cols <- c("start", "duration", "when_inf", "reached", "dir_rec")
  res_dt <- lapply(res, \(x) t(x[cols])) |>
    do.call(rbind, args = _) |>
    as.data.table()
  res_dt[
    , `:=`(start = lapply(start, \(x) unlist(x, FALSE, FALSE)),
           duration = lapply(duration, \(x) unlist(x, FALSE, FALSE)),
           reached = lapply(reached, \(x) 1L * x),
           dir_rec = lapply(dir_rec, \(x) 1L * x))][
    , `:=`(psk = psk, spr = spr, rec = rec,
           att_rate = unlist(lapply(reached, sum), FALSE, FALSE))]
  res_dt$max_inf <- unlist(lapply(res, \(x) max(x$n_inf)), FALSE, FALSE)
  return(res_dt)
}

get_evolution <- function(start, psk, spr, rec, func, n_sim = 100) {
  res <- lapply(
    seq_len(n_sim),
    function(i) {
      if (func == "base") {
        rumour_base(
          graph = g, n_iters = 1e4, inf_0 = start,
          p_skep = psk, spr_rate = spr, rec_rate = rec,
          seed = FALSE, display = FALSE
        )
      } else {
        rumour_dose(
          graph = g, n_iters = 1e4, inf_0 = start,
          p_skep = psk, spr_rate = spr, rec_rate = rec, thresh = 5,
          seed = FALSE, display = FALSE
        )
      }
    }
  )
  return(res[c("n_sus", "n_inf", "n_rec")])
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
    rates, \(x) get_results(0, 0, x[1], x[2], func = "dose", n_sim = 500),
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
    \(x) get_results(0, 0, x[1], x[2], func = "dose", n_sim = 500)
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}

fs_time - st_time

str(results)
# dput(results, "simulations/dose.txt")


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
    \(s, r) get_results(0, 0, s, r, func = "dose", n_sim = 500),
    SIMPLIFY = FALSE
  )
)
dim(results) <- NULL # convert to list
fs_time <- Sys.time()
fs_time - st_time

# write to file
# dput(results, "base.txt") # or "dose.txt"
# read back with
# res <- dget("out.txt")

#######################
# SKEPTICISM ANALYSIS #
#######################

# let's choose spr = 0.85, rec = 0.15
res <- rumour_dose(
  graph = g, n_iters = 1e4, inf_0 = 0,
  p_skep = 0, spr_rate = 0.85, rec_rate = 0.15, thresh = 5,
  seed = FALSE, display = TRUE
)

n_cores <- 3
p_skep <- seq(0.05, 0.95, 0.05)

if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  st_time <- Sys.time()
  results <- mclapply(
    p_skep, \(x) get_results(
      start = 0, psk = x, spr = 0.85, rec = 0.15,
      func = "base", n_sim = 1000
    ),
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  fs_time <- Sys.time()
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster, c("get_results", "rumour_base", "rumour_dose", "p_skep", "g")
  )
  st_time <- Sys.time()
  results <- parLapplyLB(
    cluster, p_skep, \(x) get_results(
      start = 0, psk = x, spr = 0.85, rec = 0.15,
      func = "base", n_sim = 1000
    )
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}

res_dt <- rbindlist(results)
setcolorder(
  res_dt, c(
    "psk", "spr", "rec", "start", "duration",
    "att_rate", "max_inf", "when_inf", "reached", "dir_rec"
  )
)

# fwrite(res_dt, "out/base_skeptics.csv")

# res_dt <- fread("out/base_skeptics.csv")
# chr_to_bool <- \(x) 1L * as.logical(tstrsplit(x, "|", fixed = TRUE))
# chr_to_int <- \(x) as.integer(tstrsplit(x, "|", fixed = TRUE))
# res_dt[, `:=`(reached = lapply(reached, chr_to_bool),
#               dir_rec = lapply(dir_rec, chr_to_bool),
#               when_inf = lapply(when_inf, chr_to_int))]
# fwrite(res_dt, "out/base_skeptics.csv")

# res_dt <- fread("out/base_skeptics.csv")
# chr_to_int <- \(x) as.integer(tstrsplit(x, "|", fixed = TRUE))
# cols <- c("reached", "dir_rec", "when_inf")
# res_dt[, (cols) := lapply(.SD, \(x) lapply(x, chr_to_int)), .SDcols = cols]

################
# NODE BY NODE #
################

st_node <- 1:vcount(g)
if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  st_time <- Sys.time()
  results <- mclapply(
    st_node, \(x) get_results(
      start = x, psk = 0.5, spr = 0.85, rec = 0.15,
      func = "base", n_sim = 10
    ),
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  fs_time <- Sys.time()
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster, c("get_results", "rumour_base", "rumour_dose", "st_node", "g")
  )
  st_time <- Sys.time()
  results <- parLapplyLB(
    cluster, st_node, \(x) get_results(
      start = x, psk = 0.5, spr = 0.85, rec = 0.15,
      func = "base", n_sim = 1000
    )
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}

res_dt <- rbindlist(results)
setcolorder(
  res_dt, c(
    "psk", "spr", "rec", "start", "duration",
    "att_rate", "max_inf", "when_inf", "reached", "dir_rec"
  )
)

# fwrite(res_dt, "out/node_by_node.csv")

######################
# OUTBREAK EVOLUTION #
######################

func <- "base"
start <- 0
psk <- 0.5
spr <- 0.85
rec <- 0.15
iters <- 1000

st_time <- Sys.time()
results <- get_evolution(start, psk, spr, rec, func, iters)
fs_time <- Sys.time()

fs_time - st_time
str(results)

dput(
  results,
  sprintf("simulations/%s_%i_%g_%g_%g.txt", func, start, psk, spr, rec)
)
