parse_ameras_formula <- function(formula, data, family) {
  specials <- c("dose", "strata", "offset")
  X_formula <- collect_X(formula)

  if (family == "prophaz") {
    surv <- parse_surv_term(formula)
    formula[[2]] <- quote(.response)
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data)

    return(list(
      Y = NULL,
      status = surv$status,
      entry = surv$entry,
      exit = surv$exit,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      X_formula = X_formula,
      offset = NULL,
      setnr = NULL
    ))
  } else if (family == "clogit") {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data)
    strata <- parse_strata_term(tt)

    if (is.null(strata$setnr)) {
      stop("Formula for family='clogit' must contain a strata() term")
    }

    return(list(
      Y = as.character(formula[[2]]),
      status = NULL,
      entry = NULL,
      exit = NULL,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      X_formula = X_formula,
      offset = NULL,
      setnr = strata$setnr
    ))
  } else if (family == "poisson") {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data)
    off <- parse_offset_term(tt)

    return(list(
      Y = as.character(formula[[2]]),
      status = NULL,
      entry = NULL,
      exit = NULL,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      X_formula = X_formula,
      offset = off$offset,
      setnr = NULL
    ))
  } else {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data)

    return(list(
      Y = as.character(formula[[2]]),
      status = NULL,
      entry = NULL,
      exit = NULL,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      X_formula = X_formula,
      offset = NULL,
      setnr = NULL
    ))
  }
}


parse_dose_term <- function(tt, data) {
  special_idx <- attr(tt, "specials")$dose
  if (is.null(special_idx)) {
    stop("Formula must contain a dose() term")
  }

  rhs <- attr(tt, "variables")
  dose_call <- rhs[[special_idx + 1]]
  dose_args <- as.list(dose_call[-1])

  nms <- names(dose_args)
  if (is.null(nms)) {
    nms <- rep("", length(dose_args))
  }
  named_idx <- !is.na(nms) & nzchar(nms)
  named_args <- dose_args[named_idx]
  sel_args <- dose_args[!named_idx]

  dosevars <- resolve_dose_selection(sel_args, data)
  doseRRmod <- if (!is.null(named_args$model)) {
    as.character(named_args$model)
  } else {
    "ERR"
  }
  deg <- if (!is.null(named_args$deg)) as.integer(named_args$deg) else 1
  M <- if (!is.null(named_args$modifier)) {
    parse_modifier(named_args$modifier)
  } else {
    NULL
  }

  list(
    dosevars = dosevars,
    doseRRmod = doseRRmod,
    deg = deg,
    M = M
  )
}


parse_modifier <- function(expr) {
  if (is.null(expr)) {
    return(NULL)
  }
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
  rhs <- formula[[3]]

  if (has_no_intercept(rhs)) {
    stop(
      "Removing the intercept via -1 or +0 is not currently supported. "
    )
  }

  remove_specials <- function(expr) {
    if (is.symbol(expr)) {
      return(expr)
    }
    if (is.call(expr)) {
      fn <- as.character(expr[[1]])
      if (fn %in% specials) {
        return(NULL)
      }
      if (fn %in% c("+", "-", "*", ":", "^")) {
        args <- lapply(as.list(expr)[-1], remove_specials)
        args <- Filter(Negate(is.null), args)
        if (!length(args)) {
          return(NULL)
        }
        if (length(args) == 1) {
          return(args[[1]])
        }
        return(as.call(c(list(expr[[1]]), args)))
      }
    }
    expr
  }

  cleaned_rhs <- remove_specials(rhs)
  if (is.null(cleaned_rhs)) {
    return(NULL)
  }

  # Return the cleaned RHS as a formula for later use by model.matrix
  as.formula(paste("~", deparse(cleaned_rhs, width.cutoff = 500L)))
}


has_no_intercept <- function(expr) {
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    if (fn == "-") {
      # Check for -1
      args <- as.list(expr)[-1]
      if (any(sapply(args, function(a) is.numeric(a) && a == 1))) {
        return(TRUE)
      }
    }
    if (fn == "+") {
      # Check for +0
      args <- as.list(expr)[-1]
      if (any(sapply(args, function(a) is.numeric(a) && a == 0))) {
        return(TRUE)
      }
    }
    # Recurse into sub-expressions
    return(any(sapply(as.list(expr)[-1], has_no_intercept)))
  }
  FALSE
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
      entry = NULL,
      exit = as.character(args[[1]]),
      status = as.character(args[[2]])
    )
  } else if (length(args) == 3) {
    # Surv(entry, exit, status)
    list(
      entry = as.character(args[[1]]),
      exit = as.character(args[[2]]),
      status = as.character(args[[3]])
    )
  } else {
    stop(
      "Surv() must have either 2 arguments (exit, status) or ",
      "3 arguments (entry, exit, status)"
    )
  }
}

parse_strata_term <- function(tt) {
  special_idx <- attr(tt, "specials")$strata
  if (is.null(special_idx)) {
    return(list(setnr = NULL))
  }

  rhs <- attr(tt, "variables")
  strata_call <- rhs[[special_idx + 1]]
  strata_args <- as.list(strata_call[-1])

  if (length(strata_args) != 1) {
    stop("strata() must contain exactly one variable")
  }

  list(setnr = as.character(strata_args[[1]]))
}


parse_offset_term <- function(tt) {
  special_idx <- attr(tt, "specials")$offset
  if (is.null(special_idx)) {
    return(list(offset = NULL))
  }

  rhs <- attr(tt, "variables")
  offset_call <- rhs[[special_idx + 1]]
  offset_args <- as.list(offset_call[-1])

  if (length(offset_args) != 1) {
    stop("offset() must contain exactly one variable")
  }

  list(offset = as.character(offset_args[[1]]))
}
