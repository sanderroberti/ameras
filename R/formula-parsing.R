parse_ameras_formula <- function(formula, data, family, env = parent.frame()) {
  specials <- c("dose", "strata", "offset")

  # Validate the formula structure before extracting model inputs. Survival
  # responses need parsing, whereas other families require one response column.
  validate_formula(formula)
  surv <- NULL
  if (family == "prophaz") {
    surv <- parse_surv_term(formula)
  } else {
    validate_response_term(formula)
  }
  validate_special_terms(formula, family)

  # terms() identifies supported special calls without evaluating them. Replace
  # the Surv() response with a placeholder while retaining the original RHS.
  terms_formula <- formula
  if (family == "prophaz") {
    terms_formula[[2]] <- quote(.response)
  }
  tt <- terms(terms_formula, specials = specials)

  X_formula <- collect_X(formula, env = env)

  if (family == "prophaz") {
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
    dose <- parse_dose_term(tt, data, env = env)
    strata <- parse_strata_term(tt)

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


# Validate properties shared by every supported model formula.
validate_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  if (length(formula) != 3L) {
    stop(
      "Model formula must be two-sided, with a response on the left-hand ",
      "side.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


# Non-survival families currently support only a single response data column.
validate_response_term <- function(formula) {

  response <- formula[[2]]
  if (!is.symbol(response)) {
    stop(
      "The response must be a single, untransformed data-column name.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


# Walk the formula syntax tree instead of relying only on terms(). This retains
# duplicate identical calls, which terms() may collapse into a single term.
find_formula_calls <- function(expr, name) {
  if (!is.call(expr)) {
    return(list())
  }

  current <- if (is.symbol(expr[[1]]) && as.character(expr[[1]]) == name) {
    list(expr)
  } else {
    list()
  }
  children <- unlist(
    lapply(as.list(expr)[-1], find_formula_calls, name = name),
    recursive = FALSE
  )

  c(current, children)
}


# A special is standalone when it appears directly in the additive formula
# structure, rather than inside an interaction, transformation, or function.
is_standalone_formula_term <- function(expr, name) {
  if (!is.call(expr)) {
    return(FALSE)
  }

  fn <- if (is.symbol(expr[[1]])) as.character(expr[[1]]) else ""
  if (fn == name) {
    return(TRUE)
  }
  if (fn %in% c("+", "(")) {
    return(any(vapply(
      as.list(expr)[-1],
      is_standalone_formula_term,
      logical(1),
      name = name
    )))
  }
  if (fn == "-" && length(expr) == 3L) {
    return(is_standalone_formula_term(expr[[2]], name))
  }

  FALSE
}


# Check the number, placement, and family-specific availability of special
# terms before collect_X() removes them from the ordinary covariate formula.
validate_special_terms <- function(formula, family) {
  rhs <- formula[[3]]
  dose_calls <- find_formula_calls(rhs, "dose")
  strata_calls <- find_formula_calls(rhs, "strata")
  offset_calls <- find_formula_calls(rhs, "offset")

  if (length(dose_calls) != 1L) {
    stop(
      "Formula must contain exactly one standalone dose() term.",
      call. = FALSE
    )
  }
  if (!is_standalone_formula_term(rhs, "dose")) {
    stop(
      "dose() must be a standalone additive term and cannot be used in an ",
      "interaction or transformation.",
      call. = FALSE
    )
  }

  if (family == "clogit") {
    if (length(strata_calls) != 1L) {
      stop(
        "Formula for family='clogit' must contain exactly one strata() term.",
        call. = FALSE
      )
    }
    if (!is_standalone_formula_term(rhs, "strata")) {
      stop(
        "strata() must be a standalone additive term.",
        call. = FALSE
      )
    }
  } else if (length(strata_calls)) {
    stop(
      "strata() is only supported for family='clogit'.",
      call. = FALSE
    )
  }

  if (family == "poisson") {
    if (length(offset_calls) > 1L) {
      stop(
        "Formula for family='poisson' can contain at most one offset() term.",
        call. = FALSE
      )
    }
    if (
      length(offset_calls) == 1L &&
        !is_standalone_formula_term(rhs, "offset")
    ) {
      stop(
        "offset() must be a standalone additive term.",
        call. = FALSE
      )
    }
  } else if (length(offset_calls)) {
    stop(
      "offset() is only supported for family='poisson'.",
      call. = FALSE
    )
  }

  NULL
}


parse_dose_term <- function(tt, data, env = parent.frame()) {
  special_idx <- attr(tt, "specials")$dose
  if (length(special_idx) != 1L) {
    stop(
      "dose() must be a standalone additive term in the formula.",
      call. = FALSE
    )
  }

  # The factors matrix records whether dose() participates in an interaction.
  # This complements the syntax-tree check above after terms() normalization.
  factors <- attr(tt, "factors")
  dose_terms <- which(factors[special_idx, , drop = TRUE] > 0)
  if (
    length(dose_terms) != 1L ||
      any(colSums(factors[, dose_terms, drop = FALSE] > 0) != 1L)
  ) {
    stop(
      "dose() must be a standalone additive term and cannot be used in an ",
      "interaction or transformation.",
      call. = FALSE
    )
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

  # dose() is parsed rather than evaluated, so reproduce strict argument
  # matching explicitly and reject misspellings or duplicate settings.
  allowed_args <- c("model", "deg", "modifier")
  unknown_args <- unique(nms[named_idx & !nms %in% allowed_args])
  if (length(unknown_args)) {
    stop(
      "Unknown named argument(s) in dose(): ",
      paste(unknown_args, collapse = ", "),
      ". Supported named arguments are model, deg, and modifier.",
      call. = FALSE
    )
  }

  duplicated_args <- unique(nms[named_idx & duplicated(nms)])
  if (length(duplicated_args)) {
    stop(
      "Argument(s) specified more than once in dose(): ",
      paste(duplicated_args, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!length(sel_args)) {
    stop(
      "dose() must contain at least one unnamed dose-selection expression.",
      call. = FALSE
    )
  }

  dosevars <- resolve_dose_selection(sel_args, data, env = env)
  doseRRmod <- if (!is.null(named_args$model)) {
    model <- named_args$model
    if (!(is.symbol(model) || (is.character(model) && length(model) == 1L))) {
      stop(
        "The model argument in dose() must be one of 'ERR', 'EXP', or ",
        "'LINEXP'.",
        call. = FALSE
      )
    }
    model <- as.character(model)
    if (!model %in% c("ERR", "EXP", "LINEXP")) {
      stop(
        "The model argument in dose() must be one of 'ERR', 'EXP', or ",
        "'LINEXP'.",
        call. = FALSE
      )
    }
    model
  } else {
    "ERR"
  }
  deg <- if (!is.null(named_args$deg)) {
    deg <- named_args$deg
    if (!is.numeric(deg) || length(deg) != 1L || !deg %in% c(1, 2)) {
      stop(
        "The deg argument in dose() must be the number 1 or 2.",
        call. = FALSE
      )
    }
    as.integer(deg)
  } else {
    1L
  }
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
  # Multiple unnamed expressions are treated like c(expr1, expr2, ...), while
  # a single expression is passed to tidyselect unchanged.
  sel_expr <- if (length(sel_args) == 1) {
    sel_args[[1]]
  } else {
    as.call(c(list(quote(c)), sel_args))
  }

  idx <- tidyselect::eval_select(sel_expr, data, env = env)
  if (!length(idx)) {
    stop(
      "The dose() selection must select at least one data column.",
      call. = FALSE
    )
  }
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
      "X formula terms must preserve one value per row of data. ",
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
        "could not build the X model matrix. X formula terms must ",
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
      "X model matrix has ",
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
  # Surv() is not evaluated here; store its plain column-name arguments for the
  # family-specific input checks and likelihood code used later.
  lhs <- formula[[2]]

  if (
    !is.call(lhs) ||
      !is.symbol(lhs[[1]]) ||
      as.character(lhs[[1]]) != "Surv"
  ) {
    stop(
      "Left-hand side must be a Surv() term for family='prophaz'.",
      call. = FALSE
    )
  }

  args <- as.list(lhs[-1])

  if (!length(args) %in% c(2L, 3L)) {
    stop(
      "Surv() must have either 2 arguments (exit, status) or ",
      "3 arguments (entry, exit, status).",
      call. = FALSE
    )
  }
  if (!all(vapply(args, is.symbol, logical(1)))) {
    stop(
      "Every Surv() argument must be a single, untransformed data-column ",
      "name.",
      call. = FALSE
    )
  }

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
  }
}

parse_strata_term <- function(tt) {
  # terms() supplies the unevaluated special call, from which we retain only
  # the matched-set column name.
  special_idx <- attr(tt, "specials")$strata
  if (is.null(special_idx)) {
    return(list(setnr = NULL))
  }

  rhs <- attr(tt, "variables")
  strata_call <- rhs[[special_idx + 1]]
  strata_args <- as.list(strata_call[-1])

  if (length(strata_args) != 1L || !is.symbol(strata_args[[1]])) {
    stop(
      "strata() must contain exactly one untransformed data-column name.",
      call. = FALSE
    )
  }

  list(setnr = as.character(strata_args[[1]]))
}


parse_offset_term <- function(tt) {
  # As with strata(), offsets are stored as column names rather than evaluated
  # formula expressions.
  special_idx <- attr(tt, "specials")$offset
  if (is.null(special_idx)) {
    return(list(offset = NULL))
  }

  rhs <- attr(tt, "variables")
  offset_call <- rhs[[special_idx + 1]]
  offset_args <- as.list(offset_call[-1])

  if (length(offset_args) != 1L || !is.symbol(offset_args[[1]])) {
    stop(
      "offset() must contain exactly one untransformed data-column name.",
      call. = FALSE
    )
  }

  list(offset = as.character(offset_args[[1]]))
}
