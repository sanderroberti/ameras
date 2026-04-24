snapshot_combinations <- list(
  list(family = "binomial", Y = "Y.binomial", doseRRmod = "ERR", deg = 1),
  list(family = "binomial", Y = "Y.binomial", doseRRmod = "ERR", deg = 2),
  list(family = "binomial", Y = "Y.binomial", doseRRmod = "EXP", deg = 1),
  list(family = "binomial", Y = "Y.binomial", doseRRmod = "EXP", deg = 2),
  list(family = "binomial", Y = "Y.binomial", doseRRmod = "LINEXP", deg = 2),
  list(family = "poisson", Y = "Y.poisson", doseRRmod = "ERR", deg = 1),
  list(family = "poisson", Y = "Y.poisson", doseRRmod = "EXP", deg = 2),
  list(family = "gaussian", Y = "Y.gaussian", doseRRmod = NULL, deg = 1),
  list(family = "gaussian", Y = "Y.gaussian", doseRRmod = NULL, deg = 2),
  list(family = "prophaz", Y = "status", doseRRmod = "EXP", deg = 1),
  list(family = "multinomial", Y = "Y.multinomial", doseRRmod = "EXP", deg = 1),
  list(family = "clogit", Y = "Y.clogit", doseRRmod = "EXP", deg = 1)
)

covariate_combinations <- list(
  list(X = NULL, M = NULL),
  list(X = "X1", M = NULL),
  list(X = c("X1", "X2"), M = NULL),
  list(X = NULL, M = "M1"),
  list(X = "X1", M = "M1"),
  list(X = c("X1", "X2"), M = "M1"),
  list(X = NULL, M = c("M1", "M2")),
  list(X = "X1", M = c("M1", "M2")),
  list(X = c("X1", "X2"), M = c("M1", "M2"))
)

all_methods <- c("RC", "ERC", "MCML", "FMA", "BMA")

fit_combination <- function(
  family,
  Y,
  doseRRmod = NULL,
  deg,
  X,
  M,
  methods = "RC",
  data,
  niter.BMA = 1000,
  nburnin.BMA = 200,
  nthin.BMA = 1,
  transform = NULL,
  transform.jacobian = NULL,
  ...
) {
  if (length(M) > 0) {
    Mstring <- paste(M, collapse = "+")
  } else {
    Mstring <- "NULL"
  }
  if (length(X) > 0) {
    Xstring <- paste("+", X, collapse = "+")
  } else {
    Xstring <- ""
  }

  if (family == "gaussian") {
    myform <- as.formula(paste0(
      Y,
      "~dose(V1:V10, modifier=",
      Mstring,
      ", deg=",
      deg,
      ")",
      Xstring
    ))
  } else if (family == "prophaz") {
    myform <- as.formula(paste0(
      "Surv(time,",
      Y,
      ")",
      "~dose(V1:V10, modifier=",
      Mstring,
      ", deg=",
      deg,
      ",model=",
      doseRRmod,
      ")",
      Xstring
    ))
  } else if (family == "clogit") {
    myform <- as.formula(paste0(
      Y,
      "~dose(V1:V10, modifier=",
      Mstring,
      ", deg=",
      deg,
      ",model=",
      doseRRmod,
      ")",
      Xstring,
      "+ strata(setnr)"
    ))
  } else {
    myform <- as.formula(paste0(
      Y,
      "~dose(V1:V10, modifier=",
      Mstring,
      ", deg=",
      deg,
      ",model=",
      doseRRmod,
      ")",
      Xstring
    ))
  }

  args <- c(
    list(
      formula = myform,
      data = data,
      family = family,
      methods = methods,
      nburnin.BMA = nburnin.BMA,
      niter.BMA = niter.BMA,
      nthin.BMA = nthin.BMA,
      transform = transform,
      transform.jacobian = transform.jacobian
    ),
    list(...)
  )

  suppressWarnings(do.call(ameras, args))
}
