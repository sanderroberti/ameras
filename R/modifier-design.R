# modifier_info tracks three related but distinct concepts:
# source_vars are user-facing variables from the formula, design_vars are the
# numeric columns used by likelihoods, and parameter_names/group_labels are what
# users see in coefficients and intervals.
new_modifier_info <- function(
  source_vars = NULL,
  coding = "none",
  design_vars = NULL,
  parameter_names = NULL,
  group_labels = NULL,
  levels = NULL
) {
  list(
    source_vars = source_vars,
    coding = coding,
    design_vars = design_vars,
    parameter_names = parameter_names,
    group_labels = group_labels,
    levels = levels
  )
}


modifier_is_group_coded <- function(modifier_info) {
  identical(modifier_info$coding, "group")
}


modifier_source_vars <- function(m) {
  # Missing-data handling and keep.data = FALSE reconstruction need the
  # original user variables, not generated numeric design columns.
  info <- m$modifier_info
  if (!is.null(info$source_vars)) {
    return(info$source_vars)
  }

  m$M_names
}


parse_modifier_formula <- function(expr) {
  # An intercept in the modifier formula means reference-plus-contrast coding
  # (~ M). Removing the intercept means subgroup-specific coding (~ 0 + M or
  # ~ M - 1).
  f <- as.formula(expr, env = parent.frame())
  tt <- terms(f)
  vars <- all.vars(f)
  term_labels <- attr(tt, "term.labels")
  has_intercept <- attr(tt, "intercept") == 1

  if (!length(vars)) {
    return(new_modifier_info())
  }
  if (!identical(term_labels, vars)) {
    stop(
      "Interactions and transformed terms in modifier formulas are not ",
      "currently supported. Use only additive binary variables, such as ",
      "modifier = ~ M1 + M2 or modifier = ~ 0 + M.",
      call. = FALSE
    )
  }

  coding <- if (has_intercept) "contrast" else "group"
  if (identical(coding, "group") && length(vars) != 1) {
    stop(
      "Subgroup-coded modifiers currently support exactly one binary ",
      "modifier variable.",
      call. = FALSE
    )
  }

  new_modifier_info(
    source_vars = vars,
    coding = coding,
    parameter_names = vars
  )
}


parse_legacy_modifier <- function(expr, env = parent.frame()) {
  if (is.call(expr)) {
    fn <- call_name(expr)
    if (fn %in% c("*", ":", "^")) {
      stop(
        "Interactions in modifier are not currently supported. ",
        "Please create interaction terms manually as new columns in your ",
        "data frame before calling ameras(). ",
        "For example: data$M1M2 <- data$M1 * data$M2",
        call. = FALSE
      )
    }
  }

  vars <- all.vars(expr)
  if (!length(vars)) {
    return(new_modifier_info())
  }

  lifecycle::deprecate_warn(
    when = "0.5.0",
    what = "dose(modifier = 'M1 + M2')",
    details = "Use formula modifier syntax, e.g. dose(..., modifier = ~ M1 + M2).",
    user_env = env
  )

  new_modifier_info(
    source_vars = vars,
    coding = "legacy",
    design_vars = vars,
    parameter_names = vars
  )
}


parse_modifier <- function(expr, env = parent.frame()) {
  if (is.null(expr) || identical(expr, quote(NULL))) {
    return(new_modifier_info())
  }
  if (is.call(expr) && identical(call_name(expr), "~")) {
    return(parse_modifier_formula(expr))
  }

  parse_legacy_modifier(expr, env = env)
}


encode_binary_modifier <- function(x, var) {
  # Store formula modifiers as 0/1 columns because the existing likelihoods
  # already operate on a numeric M matrix.
  if (is.factor(x)) {
    lv <- levels(x)
    if (length(lv) != 2) {
      stop(
        "ERROR: modifier variable ",
        var,
        " must have exactly two levels.",
        call. = FALSE
      )
    }
    return(list(values = as.integer(x == lv[2]), levels = lv))
  }

  if (is.logical(x)) {
    return(list(values = as.integer(x), levels = c("FALSE", "TRUE")))
  }

  if (is.numeric(x) && all(x %in% c(0, 1))) {
    return(list(values = as.numeric(x), levels = c("0", "1")))
  }

  stop(
    "ERROR: modifier variable ",
    var,
    " must be binary numeric/logical or a factor with exactly two levels.",
    call. = FALSE
  )
}


make_modifier_design_names <- function(data, n, existing = NULL) {
  # Generated columns are internal implementation details. Avoid reusing user
  # column names so a supplied data frame can be reconstructed unambiguously.
  if (!is.null(existing)) {
    return(existing)
  }

  out <- character(n)
  used <- colnames(data)
  for (i in seq_len(n)) {
    base <- paste0(".ameras_modifier_", i)
    nm <- base
    suffix <- 1L
    while (nm %in% c(used, out)) {
      suffix <- suffix + 1L
      nm <- paste0(base, "_", suffix)
    }
    out[i] <- nm
  }
  out
}


# Formula modifiers may need generated numeric columns before the existing
# likelihood code can use them. Return the augmented data, the completed
# modifier metadata, and the generated column names separately so callers do not
# have to infer those pieces from a single "design" object.
prepare_modifier_inputs <- function(data, modifier_info) {
  if (is.null(modifier_info) || !length(modifier_info$source_vars)) {
    return(list(
      data = data,
      modifier_info = new_modifier_info(),
      design_vars = NULL
    ))
  }

  if (identical(modifier_info$coding, "legacy")) {
    modifier_info$design_vars <- modifier_info$source_vars
    modifier_info$parameter_names <- modifier_info$source_vars
    return(list(
      data = data,
      modifier_info = modifier_info,
      design_vars = modifier_info$design_vars
    ))
  }

  design_vars <- make_modifier_design_names(
    data,
    length(modifier_info$source_vars),
    existing = modifier_info$design_vars
  )
  levels <- vector("list", length(modifier_info$source_vars))
  names(levels) <- modifier_info$source_vars

  # Formula modifiers keep source-variable names for reporting, but use
  # generated numeric columns for fitting.
  for (i in seq_along(modifier_info$source_vars)) {
    var <- modifier_info$source_vars[i]
    encoded <- encode_binary_modifier(data[[var]], var)
    data[[design_vars[i]]] <- encoded$values
    levels[[var]] <- encoded$levels
  }

  modifier_info$design_vars <- design_vars
  modifier_info$levels <- levels
  modifier_info$parameter_names <- modifier_info$source_vars

  if (modifier_is_group_coded(modifier_info)) {
    var <- modifier_info$source_vars[1]
    modifier_info$group_labels <- paste0(var, "=", levels[[var]])
  }

  list(data = data, modifier_info = modifier_info, design_vars = design_vars)
}


modifier_display_names <- function(data, M = NULL, modifier_info = NULL) {
  if (is.null(M)) {
    return(NULL)
  }
  if (!is.null(modifier_info$parameter_names)) {
    return(modifier_info$parameter_names)
  }

  names(data[, M, drop = FALSE])
}


modifier_dose_names <- function(
  doseRRmod,
  deg,
  data,
  M = NULL,
  modifier_info = NULL
) {
  if (is.null(M)) {
    if (identical(doseRRmod, "LINEXP")) {
      return(c("dose_linear", "dose_exponential"))
    }
    return(c("dose", "dose_squared")[seq_len(deg)])
  }

  if (modifier_is_group_coded(modifier_info)) {
    # Group-coded names are reported directly as subgroup effects. Internally
    # the likelihood still receives the equivalent reference/contrast
    # parameterization.
    labs <- modifier_info$group_labels
    if (identical(doseRRmod, "LINEXP")) {
      return(c(
        paste0("dose_linear[", labs[1], "]"),
        paste0("dose_exponential[", labs[1], "]"),
        paste0("dose_linear[", labs[2], "]"),
        paste0("dose_exponential[", labs[2], "]")
      ))
    }
    base <- c("dose", "dose_squared")[seq_len(deg)]
    return(c(
      paste0(base, "[", labs[1], "]"),
      paste0(base, "[", labs[2], "]")
    ))
  }

  m_names <- modifier_display_names(data, M, modifier_info)
  if (identical(doseRRmod, "LINEXP")) {
    return(c(
      c("dose_linear", "dose_exponential"),
      paste0("dose_linear:", m_names),
      paste0("dose_exponential:", m_names)
    ))
  }

  out <- c("dose", "dose_squared")[seq_len(deg)]
  out <- c(out, paste0("dose:", m_names))
  if (deg == 2) {
    out <- c(out, paste0("dose_squared:", m_names))
  }
  out
}


modifier_reported_to_internal_params <- function(
  params,
  family,
  X = NULL,
  M = NULL,
  deg = 1,
  modifier_info = NULL,
  Y = NULL,
  data = NULL
) {
  if (!modifier_is_group_coded(modifier_info) || is.null(M)) {
    return(params)
  }

  # Existing likelihood functions expect the legacy parameterization:
  #   first subgroup effect, then contrast = second - first.
  # Subgroup-coded fits report and optimize both subgroup effects directly, so
  # map only the internal contrast slot before evaluating the likelihood.
  intercept <- as.integer(!(family %in% c("prophaz", "clogit")))
  x_len <- length(X)
  block_size <- intercept + x_len + 2 * deg
  if (identical(family, "gaussian")) {
    block_size <- block_size + 1L
  }
  n_blocks <- 1L
  if (identical(family, "multinomial")) {
    n_blocks <- length(levels(data[, Y])) - 1L
    block_size <- 1L + x_len + 2 * deg
  }

  out <- params
  for (block in seq_len(n_blocks)) {
    offset <- (block - 1L) * block_size
    dose_pos <- offset + intercept + x_len + seq_len(deg)
    group2_pos <- dose_pos + deg
    out[group2_pos] <- params[group2_pos] - params[dose_pos]
  }

  out
}


make_modifier_loglik_transform <- function(
  transform = NULL,
  modifier_info = NULL,
  family,
  X = NULL,
  M = NULL,
  deg = 1,
  Y = NULL,
  data = NULL
) {
  if (!modifier_is_group_coded(modifier_info)) {
    return(transform)
  }

  function(params, ...) {
    # User/default transforms act on the reported parameter scale first. Only
    # after that do we map subgroup effects to the legacy contrast scale used by
    # the low-level likelihoods.
    if (!is.null(transform)) {
      params <- transform(params = params, ...)
    }
    modifier_reported_to_internal_params(
      params = params,
      family = family,
      X = X,
      M = M,
      deg = deg,
      modifier_info = modifier_info,
      Y = Y,
      data = data
    )
  }
}


err_default_transform_settings <- function(
  family,
  X = NULL,
  M = NULL,
  deg = 1,
  modifier_info = NULL,
  Y = NULL,
  data = NULL,
  lowlimit
) {
  # Return the index.t/lowlimit pair for the default ERR transformation. These
  # values are passed directly to make_transform() and make_transform.jacobian().
  # For group-coded ERR models, bounds apply to each reported subgroup effect.
  # This keeps default transformations on the same scale users see in coef().
  intercept <- as.integer(!(family %in% c("prophaz", "clogit")))
  x_len <- length(X)
  m_len <- length(M)
  block_size <- intercept + x_len + deg + m_len * deg
  if (identical(family, "gaussian")) {
    block_size <- block_size + 1L
  }
  n_blocks <- 1L
  if (identical(family, "multinomial")) {
    n_blocks <- length(levels(data[, Y])) - 1L
    block_size <- 1L + x_len + deg + m_len * deg
  }

  if (modifier_is_group_coded(modifier_info)) {
    within_block <- x_len + intercept + seq_len(2 * deg)
    lowlimit <- rep(lowlimit, 2)
  } else {
    within_block <- x_len + intercept + seq_len(deg)
  }

  index <- do.call(
    "c",
    lapply(seq_len(n_blocks), function(block) {
      (block - 1L) * block_size + within_block
    })
  )

  list(
    index.t = index,
    lowlimit = rep(lowlimit, n_blocks)
  )
}


linexp_default_transform_settings <- function(
  family,
  X = NULL,
  M = NULL,
  deg = 2,
  modifier_info = NULL,
  Y = NULL,
  data = NULL,
  lowlimit = 0
) {
  # Return the index.t/lowlimit pair for the default LINEXP transformation.
  # LINEXP only constrains the linear component. With subgroup coding, each
  # subgroup gets its own constrained linear component.
  intercept <- as.integer(!(family %in% c("prophaz", "clogit")))
  x_len <- length(X)
  m_len <- length(M)
  block_size <- intercept + x_len + deg + m_len * deg
  n_blocks <- 1L
  if (identical(family, "multinomial")) {
    n_blocks <- length(levels(data[, Y])) - 1L
    block_size <- 1L + x_len + deg + m_len * deg
  }

  within_block <- x_len + intercept + 1L
  if (modifier_is_group_coded(modifier_info)) {
    within_block <- c(within_block, within_block + deg)
    lowlimit <- rep(lowlimit, 2)
  }

  index <- do.call(
    "c",
    lapply(seq_len(n_blocks), function(block) {
      (block - 1L) * block_size + within_block
    })
  )

  list(
    index.t = index,
    lowlimit = rep(lowlimit, n_blocks)
  )
}
