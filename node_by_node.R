library(igraph)
library(ggplot2)
library(data.table)
library(parallel)
source("src/utils.R")
source("src/rumour_base.R")
source("src/rumour_dose.R")

if (Sys.info()["sysname"] == "Darwin") {
  setwd("/Users/lorenzobarbiero/Documents/GitHub/conspiracy-networks/")
}

g <- read_graph("data/graph_cph.graphml", format = "graphml")
k <- unname(degree(g))
n_nodes <- vcount(g)

get_everything <- function(psk, spr, rec, func) {
  res <- lapply(
    seq_len(n_nodes),
    function(i) {
      if (func == "base") {
        rumour_base(
          graph = g, n_iters = 1e4, inf_0 = i,
          p_skep = psk, spr_rate = spr, rec_rate = rec
        )
      } else {
        rumour_dose(
          graph = g, n_iters = 1e4, inf_0 = i,
          p_skep = psk, spr_rate = spr, rec_rate = rec, thresh = 5,
        )
      }
    }
  )
}

n_sim <- 2
func <- "dose"
n_cores <- 3
if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  st_time <- Sys.time()
  results <- mclapply(
    seq_len(n_sim),
    function(i) {
      res <- get_everything(psk = 0.1, spr = 0.85, rec = 0.15, func = func)
      message(sprintf("Simulation %d / %d", i, n_sim))
      dput(res, sprintf("out/nbn/%s-%03d.txt", func, i))
      return(res)
    },
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  fs_time <- Sys.time()
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster,
    c(
      "get_everything", "rumour_base", "rumour_dose",
      "n_sim", "g", "n_nodes", "func"
    )
  )
  st_time <- Sys.time()
  results <- parLapplyLB(
    cluster, seq_len(n_sim),
    function(i) {
      res <- get_everything(psk = 0.1, spr = 0.85, rec = 0.15, func = func)
      message(sprintf("Simulation %d / %d", i, n_sim))
      dput(res, sprintf("out/nbn/%s-%03d.txt", func, i))
      return(res)
    }
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}
