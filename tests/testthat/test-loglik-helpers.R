make_loglik_test_object <- function(
  family,
  Y = NULL,
  status = NULL,
  entry = NULL,
  exit = NULL,
  setnr = NULL,
  doseRRmod = NULL,
  deg = 1,
  offset = NULL
) {
  list(
    model = list(
      family = family,
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
    transform.args = list()
  )
}

expect_single_loglik_equal <- function(object, params, direct) {
  data("data", package = "ameras")
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
