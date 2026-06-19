## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

future_available <- requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)
run_parallel_examples <- future_available

## ----setup--------------------------------------------------------------------
library(ameras)

## ----data---------------------------------------------------------------------
data(data, package = "ameras")

## ----sequential-plan, eval = FALSE--------------------------------------------
# future::plan(future::sequential)
# 
# fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
#               data = data,
#               family = "gaussian",
#               methods = "FMA")

## ----multisession-plan, eval = FALSE------------------------------------------
# future::plan(future::multisession, workers = 2)
# 
# fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
#               data = data,
#               family = "gaussian",
#               methods = "FMA")
# 
# future::plan(future::sequential)

## ----chunk-size, eval = FALSE-------------------------------------------------
# future::plan(future::multisession, workers = 2)
# 
# fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
#               data = data,
#               family = "gaussian",
#               methods = "FMA",
#               future.chunk.size.FMA = 2)
# 
# future::plan(future::sequential)

## ----compare-plans, eval = run_parallel_examples------------------------------
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)

future::plan(future::sequential)
set.seed(2024)
sequential_elapsed <- system.time(
  sequential_fit <- suppressWarnings(
    ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
           data = data,
           family = "gaussian",
           methods = "FMA",
           MFMA = 10000)
  )
)[["elapsed"]]

future::plan(future::multisession, workers = 2)
set.seed(2024)
parallel_elapsed <- system.time(
  parallel_fit <- suppressWarnings(
    ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
           data = data,
           family = "gaussian",
           methods = "FMA",
           MFMA = 10000,
           future.chunk.size.FMA = 2)
  )
)[["elapsed"]]

future::plan(old_plan)

## ----compare-results, eval = run_parallel_examples----------------------------
coefficient_comparison <- data.frame(
  term = names(sequential_fit$FMA$coefficients),
  sequential = unname(sequential_fit$FMA$coefficients),
  parallel = unname(parallel_fit$FMA$coefficients),
  row.names = NULL
)

coefficient_comparison[-1] <- round(coefficient_comparison[-1], 4)
coefficient_comparison

data.frame(
  plan = c("sequential", "multisession"),
  elapsed_seconds = round(c(sequential_elapsed, parallel_elapsed), 2)
)

## ----reproducibility, eval = FALSE--------------------------------------------
# future::plan(future::multisession, workers = 2)
# 
# set.seed(2024)
# fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
#               data = data,
#               family = "gaussian",
#               methods = "FMA")
# 
# future::plan(future::sequential)

