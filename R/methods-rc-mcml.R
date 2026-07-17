ameras.mcml <- function(
  family,
  dosevars,
  data,
  deg,
  transform = NULL,
  transform.jacobian = NULL,
  Y = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  inpar = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = setnr,
  doseRRmod = NULL,
  modifier_info = NULL,
  loglim = 1e-30,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  ...
) {
  # if(length(CI)==0) stop("No CI method specified")
  # CI <- CI[CI %in% c("wald.orig","wald.transformed", "proflik")]
  # if(length(CI)>1) stop("Provide one type of CI for MCML: one of wald.orig, wald.transformed, and proflik")
  # if(length(CI)==0) stop("Incorrect CI method specified, should be one of wald.orig, wald.transformed, and proflik")
  # if(CI=="proflik" & !(params.profCI %in% c("dose","all"))) stop("Incorrect choice of parameters for profile likelihood CI supplied, should be either all or dose")
  # if(CI=="wald.transformed" & is.null(transform)) stop("No transformation specified, specify transformation or choose a different CI type")
  #
  # if(CI=="proflik") message("Note: computation times for profile likelihood intervals for MCML may be extensive with large datasets or complex models")

  if (family == "gaussian") {
    if (is.null(Y)) {
      stop("Y is required for family=gaussian")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 2 + length(X) + length(M) * deg + deg)
    }
  } else if (family == "binomial") {
    if (is.null(Y)) {
      stop("Y is required for family=binomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }
  } else if (family == "poisson") {
    if (is.null(Y)) {
      stop("Y is required for family=poisson")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }
  } else if (family == "clogit") {
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=clogit")
    }
    if (is.null(status)) {
      stop("status is required for family=clogit")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }
  } else if (family == "prophaz") {
    if (is.null(exit)) {
      stop("exit is required for family=prophaz")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=prophaz")
    }
    if (is.null(status)) {
      stop("status is required for family=prophaz")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }
  } else if (family == "multinomial") {
    if (is.null(Y)) {
      stop("Y is required for family=multinomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(
        0,
        (length(unique(data[, Y])) - 1) *
          (1 + length(X) + length(M) * deg + deg)
      )
    }
  }
  # Group-coded modifiers are reported and optimized on the subgroup scale.
  # modifier_info lets the low-level relative-risk helper evaluate them without
  # converting through numerically unstable reference/contrast differences.
  loglik_transform <- make_modifier_loglik_transform(
    transform = transform,
    modifier_info = modifier_info,
    family = family,
    X = X,
    M = M,
    deg = deg,
    Y = Y,
    data = data
  )

  loglik.mcml <- make_mcml_loglik_fn(
    family = family,
    dosevars = dosevars,
    data = data,
    Y = Y,
    M = M,
    X = X,
    offset = offset,
    entry = entry,
    exit = exit,
    status = status,
    setnr = setnr,
    doseRRmod = doseRRmod,
    deg = deg,
    loglim = loglim,
    transform = loglik_transform,
    modifier_info = modifier_info,
    ...
  )

  parnames <- make_method_parnames(
    family = family,
    data = data,
    Y = Y,
    X = X,
    M = M,
    deg = deg,
    doseRRmod = doseRRmod,
    modifier_info = modifier_info
  )
  t0 <- proc.time()
  fit <- fit_objective_with_hessian(
    start = inpar,
    fn = loglik.mcml,
    optim.method = optim.method,
    control = control,
    use_optimize = length(parnames) == 1,
    ...
  )

  out <- assemble_frequentist_fit_result(
    fit = fit,
    parnames = parnames,
    t0 = t0,
    transform = transform,
    transform.jacobian = transform.jacobian,
    ...
  )
  return(out)
}

ameras.rc <- function(
  family,
  dosevars,
  data,
  deg,
  ERC = FALSE,
  transform = NULL,
  transform.jacobian = NULL,
  Y = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  inpar = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = NULL,
  doseRRmod = NULL,
  modifier_info = NULL,
  loglim = 1e-30,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  ...
) {
  if (ERC & family != "poisson" & family != "prophaz") {
    Kmat <- cov(t(data[, dosevars, drop = FALSE]))
  } else {
    Kmat <- NULL
  }
  # Group-coded modifiers are reported and optimized on the subgroup scale.
  # modifier_info lets the low-level relative-risk helper evaluate them without
  # converting through numerically unstable reference/contrast differences.
  loglik_transform <- make_modifier_loglik_transform(
    transform = transform,
    modifier_info = modifier_info,
    family = family,
    X = X,
    M = M,
    deg = deg,
    Y = Y,
    data = data
  )

  t0 <- proc.time()
  if (family == "gaussian") {
    if (is.null(Y)) {
      stop("Y is required for family=gaussian")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 2 + length(X) + length(M) * deg + deg)
    }

    loglik.rc <- function(params, ...) {
      loglik.gaussian(
        params,
        D = "rcdose_ameras",
        X = X,
        Y = Y,
        M = M,
        data = data,
        deg = deg,
        ERC = ERC,
        Kmat = Kmat,
        loglim = loglim,
        transform = loglik_transform,
        modifier_info = modifier_info,
        ...
      )
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = FALSE,
      ...
    )

  } else if (family == "binomial") {
    if (is.null(Y)) {
      stop("Y is required for family=binomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }

    loglik.rc <- function(params, ...) {
      loglik.binomial(
        params,
        D = "rcdose_ameras",
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = ERC,
        Kmat = Kmat,
        loglim = loglim,
        transform = loglik_transform,
        modifier_info = modifier_info,
        ...
      )
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = FALSE,
      ...
    )

  } else if (family == "poisson") {
    if (is.null(Y)) {
      stop("Y is required for family=poisson")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }
    if (ERC) {
      dosemat_poisson <- as.matrix(data[, dosevars, drop = FALSE])
      storage.mode(dosemat_poisson) <- "double"
      Xc <- dosemat_poisson - data[, "rcdose_ameras"]
      Kmat_diag <- rowSums(Xc^2) / (ncol(dosemat_poisson) - 1)
      rm(dosemat_poisson)
      gc()
    }

    loglik.rc <- function(params, ...) {
      if (ERC) {
        loglik.poisson.erc(
          params,
          D = dosevars,
          X = X,
          Y = Y,
          offset = offset,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          loglim = loglim,
          transform = loglik_transform,
          modifier_info = modifier_info,
          Xc = Xc,
          Kmat_diag = Kmat_diag,
          ...
        )
      } else {
        loglik.poisson(
          params,
          D = "rcdose_ameras",
          X = X,
          Y = Y,
          offset = offset,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          loglim = loglim,
          transform = loglik_transform,
          modifier_info = modifier_info,
          ...
        )
      }
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = FALSE,
      ...
    )
  } else if (family %in% c("clogit")) {
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=clogit")
    }
    if (is.null(status)) {
      stop("status is required for family=clogit")
    }

    designmat <- t(model.matrix(~ as.factor(data[, setnr]) - 1))
    set_members <- lapply(sort(unique(data[, setnr])), function(s) {
      which(data[, setnr] == s) - 1L # zero-indexed for C++
    })

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }

    loglik.rc <- function(params, ...) {
      loglik.clogit(
        params,
        D = "rcdose_ameras",
        status = status,
        X = X,
        M = M,
        doseRRmod = doseRRmod,
        designmat = designmat,
        set_members = set_members,
        entry = entry,
        exit = exit,
        data = data,
        deg = deg,
        ERC = ERC,
        Kmat = Kmat,
        loglim = loglim,
        transform = loglik_transform,
        modifier_info = modifier_info,
        ...
      )
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = length(X) + length(M) * deg + deg == 1,
      ...
    )

  } else if (family == "prophaz") {
    if (is.null(exit)) {
      stop("exit is required for family=prophaz")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=prophaz")
    }
    if (is.null(status)) {
      stop("status is required for family=prophaz")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }

    if (ERC) {
      ord_exit <- order(data[[exit]])
      dosemat_ord <- as.matrix(
        data[ord_exit, dosevars, drop = FALSE]
      )
      storage.mode(dosemat_ord) <- "double"
      Xc_ord <- dosemat_ord - data[ord_exit, "rcdose_ameras"]
      Kmat_diag_ord <- rowSums(Xc_ord^2) / (ncol(dosemat_ord) - 1)
      rm(dosemat_ord)
      gc()
    }

    loglik.rc <- function(params, ...) {
      if (ERC) {
        loglik.prophaz.erc(
          params,
          D = dosevars,
          status = status,
          X = X,
          M = M,
          doseRRmod = doseRRmod,
          entry = entry,
          exit = exit,
          data = data,
          deg = deg,
          loglim = loglim,
          transform = loglik_transform,
          modifier_info = modifier_info,
          Xc_ord = Xc_ord,
          Kmat_diag_ord = Kmat_diag_ord,
          ...
        )
      } else {
        loglik.prophaz(
          params,
          D = "rcdose_ameras",
          status = status,
          X = X,
          M = M,
          doseRRmod = doseRRmod,
          entry = entry,
          exit = exit,
          data = data,
          deg = deg,
          loglim = loglim,
          transform = loglik_transform,
          modifier_info = modifier_info,
          ...
        )
      }
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = length(X) + length(M) * deg + deg == 1,
      ...
    )

  } else if (family == "multinomial") {
    if (is.null(Y)) {
      stop("Y is required for family=multinomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(
        0,
        (length(unique(data[, Y])) - 1) *
          (1 + length(X) + length(M) * deg + deg)
      )
    }

    loglik.rc <- function(params, ...) {
      loglik.multinomial(
        params,
        D = "rcdose_ameras",
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = ERC,
        Kmat = Kmat,
        loglim = loglim,
        transform = loglik_transform,
        modifier_info = modifier_info,
        ...
      )
    }
    fit <- fit_objective_with_hessian(
      start = inpar,
      fn = loglik.rc,
      optim.method = optim.method,
      control = control,
      use_optimize = FALSE,
      ...
    )

  }
  parnames <- make_method_parnames(
    family = family,
    data = data,
    Y = Y,
    X = X,
    M = M,
    deg = deg,
    doseRRmod = doseRRmod,
    modifier_info = modifier_info
  )

  out <- assemble_frequentist_fit_result(
    fit = fit,
    parnames = parnames,
    t0 = t0,
    transform = transform,
    transform.jacobian = transform.jacobian,
    extra = list(ERC = ERC),
    ...
  )

  return(out)
}
