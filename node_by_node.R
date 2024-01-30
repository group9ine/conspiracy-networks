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
k <- unname(degree(g))
n_nodes <- vcount(g)

get_everything <- function(psk, spr, rec, func) {
  lapply(
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

whoami <- "gg"
n_sim <- 9
done <- list.files("out/nbn") |>
  stringr::str_extract("[1-9]+") |>
  as.integer() |>
  max()
sims <- seq(done + 1, len = n_sim)

func <- "base"
n_cores <- 3

if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  st_time <- Sys.time()
  results <- mclapply(
    sims, function(i) {
      res <- get_everything(psk = 0.1, spr = 0.85, rec = 0.15, func = func)
      message(sprintf("Simulation %d/%d", i, n_sim))
      dput(res, sprintf("out/nbn/%s-%s-%03d.txt", whoami, func, i))
      return(res)
    },
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  fs_time <- Sys.time()
} else {
  cluster <- makeCluster(n_cores)
  clusterExport(
    cluster, c(
      "get_everything", "rumour_base", "rumour_dose",
      "sims", "n_sim", "g", "n_nodes", "func", "whoami"
    )
  )
  st_time <- Sys.time()
  results <- parLapplyLB(
    cluster, sims, function(i) {
      res <- get_everything(psk = 0.1, spr = 0.85, rec = 0.15, func = func)
      message(sprintf("Simulation %d/%d", i, n_sim))
      dput(res, sprintf("out/nbn/%s-%s-%03d.txt", whoami, func, i))
      return(res)
    }
  )
  fs_time <- Sys.time()
  stopCluster(cluster)
}
