make_loglik_test_object <- function(
  family,
  Y = NULL,
  status = NULL,
  entry = NULL,
  exit = NULL,
  setnr = NULL,
  doseRRmod = NULL,
  deg = 1,
  offset = NULL,
  dosevars = paste0("V", 1:10)
) {
  list(
    model = list(
      family = family,
      dosevars = dosevars,
      Y = Y,
      X = NULL,
      M = NULL,
      offset = offset,
      entry = entry,
      exit = exit,
      status = status,
      setnr = setnr,
      doseRRmod = doseRRmod,
      deg = deg,
      loglim = 1e-30
    ),
    transform = NULL,
    other.args = list()
  )
}

make_loglik_test_data <- function(dosevars = paste0("V", 1:10)) {
  data("data", package = "ameras")
  data$rcdose_ameras <- rowMeans(data[, dosevars, drop = FALSE])
  data
}

expect_single_loglik_equal <- function(object, params, direct) {
  data <- make_loglik_test_data(object$model$dosevars)
  helper <- ameras:::make_base_loglik_fn_single(
    object = object,
    method_fit = NULL,
    data = data,
    dose.col = "V1"
  )

  expect_equal(helper(params), direct(data), tolerance = 1e-12)
}

test_that("single-realization loglik helper forwards runtime arguments", {
  data("data", package = "ameras")
  transform_shift <- function(params, shift, ...) params + shift
  helper <- ameras:::make_single_realization_loglik(
    family = "gaussian",
    dose.col = "V1",
    data = data,
    Y = "Y.gaussian",
    transform = transform_shift
  )

  expect_equal(
    helper(c(0, 0.1, 1), shift = c(0.2, 0.3, 0.4)),
    ameras:::loglik.gaussian(
      params = c(0.2, 0.4, 1.4),
      D = "V1",
      X = NULL,
      Y = "Y.gaussian",
      M = NULL,
      data = data,
      deg = 1,
      ERC = FALSE,
      Kmat = NULL,
      loglim = 1e-30,
      transform = NULL
    ),
    tolerance = 1e-12
  )
})

test_that("single-dose loglik helper matches direct gaussian likelihood", {
  object <- make_loglik_test_object("gaussian", Y = "Y.gaussian")

  expect_single_loglik_equal(
    object,
    params = c(0, 0.1, 1),
    direct = function(data) {
      ameras:::loglik.gaussian(
        params = c(0, 0.1, 1),
        D = "V1",
        X = NULL,
        Y = "Y.gaussian",
        M = NULL,
        data = data,
        deg = 1,
        ERC = FALSE,
        Kmat = NULL,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})

test_that("base loglik helper forwards fitted-object other.args", {
  data <- make_loglik_test_data()
  transform_shift <- function(params, shift, ...) params + shift
  object <- make_loglik_test_object(
    "poisson",
    Y = "Y.poisson",
    doseRRmod = "EXP"
  )
  object$transform <- transform_shift
  object$other.args <- list(shift = c(0.2, 0.3))

  # make_base_loglik_fn() reconstructs likelihoods from stored amerasfit
  # objects. User-supplied transform arguments are stored as other.args there.
  helper <- ameras:::make_base_loglik_fn(
    object = object,
    method_fit = list(ERC = FALSE),
    data = data
  )

  expect_equal(
    helper(c(-1, 0.1), D = "V1"),
    ameras:::loglik.poisson(
      params = c(-0.8, 0.4),
      D = "V1",
      X = NULL,
      Y = "Y.poisson",
      M = NULL,
      offset = NULL,
      doseRRmod = "EXP",
      data = data,
      deg = 1,
      loglim = 1e-30,
      transform = NULL
    ),
    tolerance = 1e-12
  )
})

test_that("base loglik helper forwards other.args for Poisson ERC", {
  data <- make_loglik_test_data()
  transform_shift <- function(params, shift, ...) params + shift
  object <- make_loglik_test_object(
    "poisson",
    Y = "Y.poisson",
    doseRRmod = "EXP"
  )
  object$transform <- transform_shift
  object$other.args <- list(shift = c(0.2, 0.3))

  helper <- ameras:::make_base_loglik_fn(
    object = object,
    method_fit = list(ERC = TRUE),
    data = data
  )

  # Poisson ERC uses precomputed centered dose-realization residuals rather
  # than the generic Kmat approximation used by several other families.
  dosemat <- as.matrix(data[, object$model$dosevars, drop = FALSE])
  storage.mode(dosemat) <- "double"
  Xc <- dosemat - data[, "rcdose_ameras"]
  Kmat_diag <- rowSums(Xc^2) / (ncol(dosemat) - 1)

  expect_equal(
    helper(c(-1, 0.1), D = "rcdose_ameras"),
    ameras:::loglik.poisson.erc(
      params = c(-0.8, 0.4),
      D = "rcdose_ameras",
      X = NULL,
      Y = "Y.poisson",
      M = NULL,
      offset = NULL,
      doseRRmod = "EXP",
      data = data,
      deg = 1,
      loglim = 1e-30,
      transform = NULL,
      Xc = Xc,
      Kmat_diag = Kmat_diag
    ),
    tolerance = 1e-12
  )
})

test_that("base loglik helper forwards other.args for proportional hazards ERC", {
  data <- make_loglik_test_data()
  transform_shift <- function(params, shift, ...) params + shift
  object <- make_loglik_test_object(
    "prophaz",
    status = "status",
    exit = "time",
    doseRRmod = "EXP"
  )
  object$transform <- transform_shift
  object$other.args <- list(shift = 0.2)

  helper <- ameras:::make_base_loglik_fn(
    object = object,
    method_fit = list(ERC = TRUE),
    data = data
  )

  # Proportional hazards ERC performs the same residual-variance calculation as
  # Poisson ERC, but on exit-time order because the risk sets are ordered.
  ord_exit <- order(data[["time"]])
  dosemat_ord <- as.matrix(data[ord_exit, object$model$dosevars, drop = FALSE])
  storage.mode(dosemat_ord) <- "double"
  Xc_ord <- dosemat_ord - data[ord_exit, "rcdose_ameras"]
  Kmat_diag_ord <- rowSums(Xc_ord^2) / (ncol(dosemat_ord) - 1)

  expect_equal(
    helper(0.1, D = "rcdose_ameras"),
    ameras:::loglik.prophaz.erc(
      params = 0.3,
      D = "rcdose_ameras",
      X = NULL,
      status = "status",
      entry = NULL,
      exit = "time",
      M = NULL,
      doseRRmod = "EXP",
      data = data,
      deg = 1,
      loglim = 1e-30,
      transform = NULL,
      Xc_ord = Xc_ord,
      Kmat_diag_ord = Kmat_diag_ord
    ),
    tolerance = 1e-12
  )
})

test_that("base loglik helper builds clogit structures from supplied data", {
  data <- make_loglik_test_data()
  object <- make_loglik_test_object(
    "clogit",
    status = "Y.clogit",
    setnr = "setnr",
    doseRRmod = "EXP"
  )

  # keep.data = FALSE users can supply data at CI time. This guards against
  # accidentally rebuilding clogit risk-set structures from object$model$data.
  object$model$data <- NULL
  helper <- ameras:::make_base_loglik_fn(
    object = object,
    method_fit = list(ERC = TRUE),
    data = data
  )

  designmat <- t(model.matrix(~ as.factor(data[, "setnr"]) - 1))
  set_members <- lapply(sort(unique(data[, "setnr"])), function(s) {
    which(data[, "setnr"] == s) - 1L
  })
  Kmat <- cov(t(data[, object$model$dosevars, drop = FALSE]))

  expect_equal(
    helper(0.1, D = "rcdose_ameras"),
    ameras:::loglik.clogit(
      params = 0.1,
      D = "rcdose_ameras",
      X = NULL,
      status = "Y.clogit",
      M = NULL,
      doseRRmod = "EXP",
      designmat = designmat,
      set_members = set_members,
      data = data,
      deg = 1,
      ERC = TRUE,
      Kmat = Kmat,
      loglim = 1e-30,
      transform = NULL
    ),
    tolerance = 1e-12
  )
})

test_that("single-dose loglik helper matches direct binomial likelihood", {
  object <- make_loglik_test_object(
    "binomial",
    Y = "Y.binomial",
    doseRRmod = "EXP"
  )

  expect_single_loglik_equal(
    object,
    params = c(-1, 0.1),
    direct = function(data) {
      ameras:::loglik.binomial(
        params = c(-1, 0.1),
        D = "V1",
        X = NULL,
        Y = "Y.binomial",
        M = NULL,
        doseRRmod = "EXP",
        data = data,
        deg = 1,
        ERC = FALSE,
        Kmat = NULL,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})

test_that("single-dose loglik helper matches direct poisson likelihood", {
  object <- make_loglik_test_object(
    "poisson",
    Y = "Y.poisson",
    doseRRmod = "EXP"
  )

  expect_single_loglik_equal(
    object,
    params = c(-1, 0.1),
    direct = function(data) {
      ameras:::loglik.poisson(
        params = c(-1, 0.1),
        D = "V1",
        X = NULL,
        Y = "Y.poisson",
        M = NULL,
        offset = NULL,
        doseRRmod = "EXP",
        data = data,
        deg = 1,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})

test_that("single-dose loglik helper matches direct multinomial likelihood", {
  object <- make_loglik_test_object(
    "multinomial",
    Y = "Y.multinomial",
    doseRRmod = "EXP"
  )
  params <- c(-0.5, 0.1, 0.2, 0.05)

  expect_single_loglik_equal(
    object,
    params = params,
    direct = function(data) {
      ameras:::loglik.multinomial(
        params = params,
        D = "V1",
        X = NULL,
        Y = "Y.multinomial",
        M = NULL,
        doseRRmod = "EXP",
        data = data,
        deg = 1,
        ERC = FALSE,
        Kmat = NULL,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})

test_that("multinomial likelihood keeps single-X coefficient blocks as matrices", {
  data <- make_loglik_test_data()

  # With one baseline covariate and more than two outcome levels, the covariate
  # coefficient block is one row by multiple columns. This guards against that
  # one-row matrix being simplified to a vector before matrix multiplication.
  expect_true(is.finite(
    ameras:::loglik.multinomial(
      params = c(-0.5, 0.2, 0.1, 0.1, -0.2, 0.05),
      D = "V1",
      X = "X1",
      Y = "Y.multinomial",
      M = NULL,
      doseRRmod = "EXP",
      data = data,
      deg = 1,
      ERC = FALSE,
      Kmat = NULL,
      loglim = 1e-30,
      transform = NULL
    )
  ))
})

test_that("single-dose loglik helper matches direct clogit likelihood", {
  object <- make_loglik_test_object(
    "clogit",
    status = "Y.clogit",
    setnr = "setnr",
    doseRRmod = "EXP"
  )

  expect_single_loglik_equal(
    object,
    params = 0.1,
    direct = function(data) {
      designmat <- t(model.matrix(~ as.factor(data[, "setnr"]) - 1))
      set_members <- lapply(sort(unique(data[, "setnr"])), function(s) {
        which(data[, "setnr"] == s) - 1L
      })

      ameras:::loglik.clogit(
        params = 0.1,
        D = "V1",
        X = NULL,
        status = "Y.clogit",
        M = NULL,
        doseRRmod = "EXP",
        designmat = designmat,
        set_members = set_members,
        data = data,
        deg = 1,
        ERC = FALSE,
        Kmat = NULL,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})

test_that("single-dose loglik helper matches direct prophaz likelihood", {
  object <- make_loglik_test_object(
    "prophaz",
    status = "status",
    exit = "time",
    doseRRmod = "EXP"
  )

  expect_single_loglik_equal(
    object,
    params = 0.1,
    direct = function(data) {
      ameras:::loglik.prophaz(
        params = 0.1,
        D = "V1",
        X = NULL,
        status = "status",
        entry = NULL,
        exit = "time",
        M = NULL,
        doseRRmod = "EXP",
        data = data,
        deg = 1,
        loglim = 1e-30,
        transform = NULL
      )
    }
  )
})
