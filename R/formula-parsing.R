parse_ameras_formula <- function(formula, data, family, env = parent.frame()) {
  specials <- c("dose", "strata", "offset")
  X_formula <- collect_X(formula, env = env)

  if (family == "prophaz") {
    surv <- parse_surv_term(formula)
    formula[[2]] <- quote(.response)
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data, env = env)

    return(list(
      Y = NULL,
      status = surv$status,
      entry = surv$entry,
      exit = surv$exit,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      modifier_info = dose$modifier_info,
      X_formula = X_formula,
      offset = NULL,
      setnr = NULL
    ))
  } else if (family == "clogit") {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data, env = env)
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
      modifier_info = dose$modifier_info,
      X_formula = X_formula,
      offset = NULL,
      setnr = strata$setnr
    ))
  } else if (family == "poisson") {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data, env = env)
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
      modifier_info = dose$modifier_info,
      X_formula = X_formula,
      offset = off$offset,
      setnr = NULL
    ))
  } else {
    tt <- terms(formula, specials = specials)
    dose <- parse_dose_term(tt, data, env = env)

    return(list(
      Y = as.character(formula[[2]]),
      status = NULL,
      entry = NULL,
      exit = NULL,
      dosevars = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg = dose$deg,
      M = dose$M,
      modifier_info = dose$modifier_info,
      X_formula = X_formula,
      offset = NULL,
      setnr = NULL
    ))
  }
}


parse_dose_term <- function(tt, data, env = parent.frame()) {
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

  dosevars <- resolve_dose_selection(sel_args, data, env = env)
  doseRRmod <- if (!is.null(named_args$model)) {
    as.character(named_args$model)
  } else {
    "ERR"
  }
  deg <- if (!is.null(named_args$deg)) as.integer(named_args$deg) else 1
  modifier_info <- if (!is.null(named_args$modifier)) {
    parse_modifier(named_args$modifier, env = env)
  } else {
    new_modifier_info()
  }

  list(
    dosevars = dosevars,
    doseRRmod = doseRRmod,
    deg = deg,
    M = modifier_info$source_vars,
    modifier_info = modifier_info
  )
}


resolve_dose_selection <- function(sel_args, data, env = parent.frame()) {
  sel_expr <- if (length(sel_args) == 1) {
    sel_args[[1]]
  } else {
    as.call(c(list(quote(c)), sel_args))
  }

  idx <- tidyselect::eval_select(sel_expr, data, env = env)
  colnames(data)[idx]
}


collect_X <- function(formula, env = parent.frame()) {
  specials <- c("dose", "strata", "offset", "Surv")
  rhs <- formula[[3]]
  formula_env <- environment(formula) %||% env

  remove_specials <- function(expr) {
    if (is.symbol(expr)) {
      return(expr)
    }
    if (is.call(expr)) {
      fn <- call_name(expr)
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
  if (has_no_intercept(cleaned_rhs)) {
    stop(
      "Removing the intercept via -1 or +0 is not currently supported. "
    )
  }

  # Return the cleaned RHS as a formula for later use by model.matrix.
  # Constructing the formula from a call preserves namespace-qualified terms
  # such as splines::ns() without round-tripping through text.
  as.formula(as.call(list(as.name("~"), cleaned_rhs)), env = formula_env)
}


call_name <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }

  # For calls like splines::ns(...), as.character(expr[[1]]) has length > 1.
  # The first element is the operator name ("::"), which is all we need when
  # deciding whether the call is a formula operator or a special.
  as.character(expr[[1]])[[1]]
}


has_no_intercept <- function(expr) {
  if (is.call(expr)) {
    fn <- call_name(expr)
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


build_X_design <- function(X_formula, data, X_design_info = NULL) {
  if (is.null(X_formula)) {
    return(list(matrix = NULL, X_design_info = NULL))
  }

  if (is.null(X_design_info)) {
    mf <- model.frame(X_formula, data = data, na.action = na.pass)
    check_X_model_frame_rows(mf, data)
    terms_obj <- terms(mf)
    X_matrix <- build_X_model_matrix(terms_obj, data = mf)
    environment(terms_obj) <- make_X_design_environment(terms_obj)
    X_design_info <- list(
      terms = terms_obj,
      contrasts = attr(X_matrix, "contrasts"),
      xlevels = .getXlevels(terms_obj, mf)
    )
  } else {
    mf <- model.frame(
      X_design_info$terms,
      data = data,
      na.action = na.pass,
      xlev = X_design_info$xlevels
    )
    check_X_model_frame_rows(mf, data)
    X_matrix <- build_X_model_matrix(
      X_design_info$terms,
      data = mf,
      contrasts.arg = X_design_info$contrasts
    )
  }

  check_X_matrix_rows(X_matrix, data)

  list(
    matrix = X_matrix[, -1, drop = FALSE],
    X_design_info = X_design_info
  )
}


check_X_model_frame_rows <- function(mf, data) {
  expected_rows <- nrow(data)
  term_rows <- vapply(mf, NROW, integer(1))
  bad_terms <- names(term_rows)[term_rows != expected_rows]

  if (length(bad_terms)) {
    stop(
      "ERROR: X formula terms must preserve one value per row of data. ",
      "The following evaluated term(s) returned the wrong number of rows: ",
      paste(bad_terms, collapse = ", "),
      ". This can happen when a term drops missing values internally; ",
      "remove or impute missing values before fitting.",
      call. = FALSE
    )
  }

  NULL
}


build_X_model_matrix <- function(terms_obj, data, contrasts.arg = NULL) {
  tryCatch(
    model.matrix(
      terms_obj,
      data = data,
      contrasts.arg = contrasts.arg
    ),
    error = function(e) {
      stop(
        "ERROR: could not build the X model matrix. X formula terms must ",
        "preserve one value per row of data; this can happen when a term ",
        "drops missing values internally. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}


check_X_matrix_rows <- function(X_matrix, data) {
  if (nrow(X_matrix) != nrow(data)) {
    stop(
      "ERROR: X model matrix has ",
      nrow(X_matrix),
      " rows, but data has ",
      nrow(data),
      " rows. This can happen when X terms drop rows because of missing ",
      "values; remove or impute missing values before fitting.",
      call. = FALSE
    )
  }

  NULL
}


make_X_design_environment <- function(terms_obj) {
  source_env <- environment(terms_obj)
  env <- new.env(parent = baseenv())
  fun_names <- unique(collect_formula_function_names(attr(terms_obj, "predvars")))

  for (fun_name in fun_names) {
    if (exists(fun_name, envir = source_env, mode = "function", inherits = TRUE)) {
      assign(
        fun_name,
        get(fun_name, envir = source_env, mode = "function", inherits = TRUE),
        envir = env
      )
    }
  }

  env
}


collect_formula_function_names <- function(expr) {
  if (!is.call(expr)) {
    return(character())
  }

  fn <- expr[[1]]
  out <- character()

  # Namespace-qualified calls, e.g. splines::ns(), are self-contained and do
  # not need a function binding in the terms environment. Unqualified calls
  # such as ns() or bs() do, otherwise keep.data = FALSE reconstruction would
  # depend on the caller's search path.
  if (is.symbol(fn)) {
    fn_name <- as.character(fn)
    formula_operators <- c(
      "(", "{", "[", "[[", "$", "@", "::", ":::",
      "+", "-", "*", "/", "^", ":", "~", "=", "<-", "list"
    )
    if (!fn_name %in% formula_operators) {
      out <- fn_name
    }
  }

  unique(c(
    out,
    unlist(
      lapply(as.list(expr)[-1], collect_formula_function_names),
      use.names = FALSE
    )
  ))
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
