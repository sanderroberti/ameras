make_base_loglik_fn <- function(object, method_fit, data) {
  m <- object$model
  ERC <- method_fit$ERC
  if (is.null(ERC)) {
    ERC <- FALSE
  }
  Kmat <- if (ERC && m$family != "poisson") {
    cov(t(data[, m$dosevars, drop = FALSE]))
  } else {
    NULL
  }

  if (m$family == "gaussian") {
    function(params, D) {
      do.call(
        loglik.gaussian,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "poisson") {
    if (ERC) {
      function(params, D) {
        do.call(
          loglik.poisson.erc,
          c(
            list(
              params = params,
              D = m$dosevars,
              X = m$X,
              Y = m$Y,
              M = m$M,
              offset = m$offset,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform
            ),
            object$other.args
          )
        )
      }
    } else {
      function(params, D) {
        do.call(
          loglik.poisson,
          c(
            list(
              params = params,
              D = D,
              X = m$X,
              Y = m$Y,
              M = m$M,
              offset = m$offset,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform
            ),
            object$other.args
          )
        )
      }
    }
  } else if (m$family == "binomial") {
    function(params, D) {
      do.call(
        loglik.binomial,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "prophaz") {
    function(params, D) {
      do.call(
        loglik.prophaz,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            status = m$status,
            entry = m$entry,
            exit = m$exit,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "clogit") {
    designmat <- t(model.matrix(~ as.factor(m$data[, m$setnr]) - 1))
    function(params, D) {
      do.call(
        loglik.clogit,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            status = m$status,
            designmat = designmat,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "multinomial") {
    function(params, D) {
      do.call(
        loglik.multinomial,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else {
    stop("Unknown family: ", m$family)
  }
}


make_loglik_fn <- function(object, method_name, method_fit, data) {
  m <- object$model
  base_fn <- make_base_loglik_fn(object, method_fit, data)

  if (method_name == "MCML") {
    function(params) {
      logliks <- base_fn(params, D = m$dosevars)
      shifted <- pmax(pmin(-logliks - max(-logliks), 7e1), -7e1)
      -log(mean(exp(shifted))) - max(-logliks)
    }
  } else if (method_name %in% c("RC", "ERC")) {
    # Poisson ERC uses dosevars directly inside loglik.poisson.erc
    # so D argument is ignored in that case, but we pass it for consistency
    function(params) {
      base_fn(params, D = "rcdose_ameras")
    }
  } else {
    stop("make_loglik_fn not applicable for method: ", method_name)
  }
}


resolve_data <- function(object, data = NULL) {
  if (!is.null(object$model$data)) {
    return(object$model$data)
  }

  if (is.null(data)) {
    stop(
      "Data not stored on object (keep.data=FALSE). ",
      "Please supply data argument to confint()."
    )
  }

  check_df(data)

  m <- object$model

  # Check dose columns separately for a clearer error message
  missing_dose <- setdiff(m$dosevars, colnames(data))
  if (length(missing_dose)) {
    stop(
      "Supplied data is missing dose columns present during fitting: ",
      paste(missing_dose, collapse = ", ")
    )
  }

  # Check remaining required variables excluding X and dose
  # X is checked separately via X_formula, dose already checked above
  needed <- setdiff(required_vars(m), c(m$dosevars, m$X))
  missing_vars <- setdiff(needed, colnames(data))
  if (length(missing_vars)) {
    stop(
      "Supplied data is missing columns present during fitting: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Check raw X variables referenced in the X formula
  if (!is.null(m$X_formula)) {
    X_raw_vars <- all.vars(m$X_formula)
    missing_X <- setdiff(X_raw_vars, colnames(data))
    if (length(missing_X)) {
      stop(
        "Supplied data is missing X variables present during fitting: ",
        paste(missing_X, collapse = ", ")
      )
    }
  }

  # Check dose columns are numeric and finite
  for (v in m$dosevars) {
    check_num_vec(data[, v, drop = TRUE], nm = paste0("dosevars:", v))
  }

  if (nrow(data) != object$num.rows) {
    warning(
      "Supplied data has ",
      nrow(data),
      " rows but original data had ",
      object$num.rows,
      " rows."
    )
  }

  # Re-expand X formula on supplied data
  if (!is.null(m$X_formula)) {
    X_matrix <- model.matrix(m$X_formula, data = data)[, -1, drop = FALSE]
    X_colnames <- colnames(X_matrix)
    new_cols <- setdiff(X_colnames, colnames(data))
    if (length(new_cols)) {
      data[, new_cols] <- X_matrix[, new_cols, drop = FALSE]
    }
  }

  data$rcdose_ameras <- rowMeans(data[, m$dosevars, drop = FALSE])

  data
}


check_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}


`%||%` <- function(x, y) if (is.null(x)) y else x
