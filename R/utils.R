
make_base_loglik_fn <- function(object, method_fit, data) {
  
  m    <- object$model
  ERC  <- method_fit$ERC
  if(is.null(ERC)) ERC <- FALSE
  Kmat <- if (ERC && m$family != "poisson") {
    cov(t(data[, m$dosevars, drop=FALSE]))
  } else {
    NULL
  }
  
  if (m$family == "gaussian") {
    function(params, D) {
      do.call(loglik.gaussian,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "poisson") {
    if (ERC) {
      function(params, D) {
        do.call(loglik.poisson.erc,
                c(list(params    = params,
                       D         = m$dosevars,
                       X         = m$X,
                       Y         = m$Y,
                       M         = m$M,
                       offset    = m$offset,
                       doseRRmod = m$doseRRmod,
                       data      = data,
                       deg       = m$deg,
                       loglim    = m$loglim,
                       transform = object$transform),
                  object$other.args))
      }
    } else {
      function(params, D) {
        do.call(loglik.poisson,
                c(list(params    = params,
                       D         = D,
                       X         = m$X,
                       Y         = m$Y,
                       M         = m$M,
                       offset    = m$offset,
                       doseRRmod = m$doseRRmod,
                       data      = data,
                       deg       = m$deg,
                       loglim    = m$loglim,
                       transform = object$transform),
                  object$other.args))
      }
    }
  } else if (m$family == "binomial") {
    function(params, D) {
      do.call(loglik.binomial,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "prophaz") {
    function(params, D) {
      do.call(loglik.prophaz,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     status    = m$status,
                     entry     = m$entry,
                     exit      = m$exit,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "clogit") {
    designmat <- t(model.matrix(~as.factor(m$data[, m$setnr]) - 1))
    function(params, D) {
      do.call(loglik.clogit,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     status    = m$status,
                     designmat = designmat,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "multinomial") {
    function(params, D) {
      do.call(loglik.multinomial,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else {
    stop("Unknown family: ", m$family)
  }
}


make_loglik_fn <- function(object, method_name, method_fit, data) {
  
  m       <- object$model
  base_fn <- make_base_loglik_fn(object, method_fit, data)
  
  if (method_name == "MCML") {
    function(params) {
      logliks <- base_fn(params, D=m$dosevars)
      shifted <- pmax(pmin(-logliks - max(-logliks), 7e1), -7e1)
      -log(mean(exp(shifted))) - max(-logliks)
    }
  } else if (method_name %in% c("RC", "ERC")) {
    # Poisson ERC uses dosevars directly inside loglik.poisson.erc
    # so D argument is ignored in that case, but we pass it for consistency
    function(params) {
      base_fn(params, D="rcdose_ameras")
    }
  } else {
    stop("make_loglik_fn not applicable for method: ", method_name)
  }
}


resolve_data <- function(object, data=NULL) {
  
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
      paste(missing_dose, collapse=", ")
    )
  }
  
  # Check remaining required variables excluding X and dose
  # X is checked separately via X_formula, dose already checked above
  needed       <- setdiff(required_vars(m), c(m$dosevars, m$X))
  missing_vars <- setdiff(needed, colnames(data))
  if (length(missing_vars)) {
    stop(
      "Supplied data is missing columns present during fitting: ",
      paste(missing_vars, collapse=", ")
    )
  }
  
  # Check raw X variables referenced in the X formula
  if (!is.null(m$X_formula)) {
    X_raw_vars <- all.vars(m$X_formula)
    missing_X  <- setdiff(X_raw_vars, colnames(data))
    if (length(missing_X)) {
      stop(
        "Supplied data is missing X variables present during fitting: ",
        paste(missing_X, collapse=", ")
      )
    }
  }
  
  # Check dose columns are numeric and finite
  for (v in m$dosevars) {
    check_num_vec(data[, v, drop=TRUE], nm=paste0("dosevars:", v))
  }
  
  if (nrow(data) != object$num.rows) {
    warning(
      "Supplied data has ", nrow(data), " rows but original data had ",
      object$num.rows, " rows."
    )
  }
  
  # Re-expand X formula on supplied data
  if (!is.null(m$X_formula)) {
    X_matrix   <- model.matrix(m$X_formula, data=data)[, -1, drop=FALSE]
    X_colnames <- colnames(X_matrix)
    new_cols   <- setdiff(X_colnames, colnames(data))
    if (length(new_cols)) {
      data[, new_cols] <- X_matrix[, new_cols, drop=FALSE]
    }
  }
  
  data$rcdose_ameras <- rowMeans(data[, m$dosevars, drop=FALSE])
  
  data
}






getCharVecStr <- function(x, sep=",") {
  
  ret <- paste0("'", x, "'")
  ret <- paste0(ret, collapse=sep)
  ret
  
}


# Function to check that an object is a string
isString <- function(obj) {
  
  if ((length(obj) == 1) && is.character(obj)) {
    ret <- TRUE
  } else {
    ret <- FALSE
  }
  
  ret
  
} # END: isString

check_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}

getVarNumbers <- function(vars, data) {
  
  if (!length(vars)) return(NULL)
  if (is.numeric(vars)) return(vars)
  cx  <- colnames(data)
  ret <- match(vars, cx)
  ret
}




parse_surv_term <- function(formula) {
  
  lhs <- formula[[2]]
  
  if (!is.call(lhs) || as.character(lhs[[1]]) != "Surv") {
    stop("Left hand side must be a Surv() term for family='prophaz'")
  }
  
  args <- as.list(lhs[-1])
  
  if (length(args) == 2) {
    # Surv(exit, status)
    list(
      entry  = NULL,
      exit   = as.character(args[[1]]),
      status = as.character(args[[2]])
    )
  } else if (length(args) == 3) {
    # Surv(entry, exit, status)
    list(
      entry  = as.character(args[[1]]),
      exit   = as.character(args[[2]]),
      status = as.character(args[[3]])
    )
  } else {
    stop("Surv() must have either 2 arguments (exit, status) or ",
         "3 arguments (entry, exit, status)")
  }
}

parse_strata_term <- function(tt) {
  
  special_idx <- attr(tt, "specials")$strata
  if (is.null(special_idx)) return(list(setnr=NULL))
  
  rhs         <- attr(tt, "variables")
  strata_call <- rhs[[special_idx + 1]]
  strata_args <- as.list(strata_call[-1])
  
  if (length(strata_args) != 1) {
    stop("strata() must contain exactly one variable")
  }
  
  list(setnr=as.character(strata_args[[1]]))
}


parse_offset_term <- function(tt) {
  
  special_idx <- attr(tt, "specials")$offset
  if (is.null(special_idx)) return(list(offset=NULL))
  
  rhs         <- attr(tt, "variables")
  offset_call <- rhs[[special_idx + 1]]
  offset_args <- as.list(offset_call[-1])
  
  if (length(offset_args) != 1) {
    stop("offset() must contain exactly one variable")
  }
  
  list(offset=as.character(offset_args[[1]]))
}



`%||%` <- function(x, y) if (is.null(x)) y else x


parse_dose_term <- function(tt, data) {
  
  special_idx <- attr(tt, "specials")$dose
  if (is.null(special_idx)) stop("Formula must contain a dose() term")
  
  rhs       <- attr(tt, "variables")
  dose_call <- rhs[[special_idx + 1]]
  dose_args <- as.list(dose_call[-1])
  
  nms <- names(dose_args)
  if (is.null(nms)) nms <- rep("", length(dose_args))
  named_idx  <- !is.na(nms) & nzchar(nms)
  named_args <- dose_args[named_idx]
  sel_args   <- dose_args[!named_idx]
  
  dosevars  <- resolve_dose_selection(sel_args, data)
  doseRRmod <- if (!is.null(named_args$model))    as.character(named_args$model)       else "ERR"
  deg       <- if (!is.null(named_args$deg))      as.integer(named_args$deg)           else 1
  M         <- if (!is.null(named_args$modifier)) parse_modifier(named_args$modifier)  else NULL
  
  list(
    dosevars  = dosevars,
    doseRRmod = doseRRmod,
    deg       = deg,
    M         = M
  )
}


parse_modifier <- function(expr) {
  
  if (is.null(expr)) return(NULL)
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    if (fn %in% c("*", ":", "^")) {
      stop(
        "Interactions in modifier are not currently supported. ",
        "Please create interaction terms manually as new columns in your ",
        "data frame before calling ameras(). ",
        "For example: data$M1M2 <- data$M1 * data$M2"
      )
    }
  }
  all.vars(expr)
}


resolve_dose_selection <- function(sel_args, data) {
  
  sel_expr <- if (length(sel_args) == 1) {
    sel_args[[1]]
  } else {
    as.call(c(list(quote(c)), sel_args))
  }
  
  idx <- tidyselect::eval_select(sel_expr, data)
  colnames(data)[idx]
}


collect_X <- function(formula) {
  
  specials <- c("dose", "strata", "offset", "Surv")
  rhs      <- formula[[3]]
  
  remove_specials <- function(expr) {
    if (is.symbol(expr)) return(expr)
    if (is.call(expr)) {
      fn <- as.character(expr[[1]])
      if (fn %in% specials) return(NULL)
      if (fn %in% c("+", "-", "*", ":", "^")) {
        args <- lapply(as.list(expr)[-1], remove_specials)
        args <- Filter(Negate(is.null), args)
        if (!length(args)) return(NULL)
        if (length(args) == 1) return(args[[1]])
        return(as.call(c(list(expr[[1]]), args)))
      }
    }
    expr
  }
  
  cleaned_rhs <- remove_specials(rhs)
  if (is.null(cleaned_rhs)) return(NULL)
  
  # Return the cleaned RHS as a formula for later use by model.matrix
  as.formula(paste("~", deparse(cleaned_rhs, width.cutoff=500L)))
}