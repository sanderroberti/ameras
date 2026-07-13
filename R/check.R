check_df <- function(x, nm = "data") {
  if (!is.data.frame(x)) {
    stop(paste0("ERROR: ", nm, " must be a data frame"))
  }
  if (!nrow(x)) {
    stop(paste0("ERROR: ", nm, " has no rows"))
  }
  if (!ncol(x)) {
    stop(paste0("ERROR: ", nm, " has no columns"))
  }
  NULL
}

check_reserved_names <- function(data, reserved = "rcdose_ameras") {
  used <- intersect(reserved, colnames(data))
  if (length(used)) {
    stop(paste0(
      "ERROR: data contains reserved ameras column name(s): ",
      getCharVecStr(used),
      ". Please rename these columns before fitting."
    ))
  }

  NULL
}

check_family <- function(x, nm = "family") {
  valid <- c(
    "gaussian",
    "binomial",
    "poisson",
    "prophaz",
    "multinomial",
    "clogit"
  )
  check_char_vec(x, nm, valid = valid, def = NULL, len = 1)
}

check_doseRRmod <- function(x, nm = "doseRRmod") {
  valid <- c("ERR", "EXP", "LINEXP")
  check_char_vec(x, nm, valid = valid, def = NULL, len = 1)
}

check_Y <- function(v, data, family) {
  nm <- "Y"
  check_vars(data, v, nm, minlen = 1, maxlen = 0)
  vec <- data[, v, drop = TRUE]
  if (family != "multinomial") {
    binary <- nonneg <- integer <- 0
    if (family %in% c("binomial", "prophaz", "clogit")) {
      binary <- 1
    } else if (family == "poisson") {
      nonneg <- 1
      integer <- 1
    }
    check_num_vec(vec, nm, binary = binary, nonneg = nonneg, integer = integer)
    NULL
  } else {
    # multinomial
    check_factor_vec(vec, nm)
    NULL
  }
}

check_D <- function(vars, data, methods) {
  nm <- "dosevars"
  check_vars(data, vars, nm, minlen = 1, maxlen = 0)
  for (v in vars) {
    vec <- data[, v, drop = TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2)
  }
  if (length(vars) == 1 & any(c("ERC", "MCML", "FMA", "BMA") %in% methods)) {
    stop(
      "Multiple exposure realizations required for ERC, MCML, FMA, and BMA. With a single exposure vector, use RC"
    )
  }
  NULL
}

check_included_realizations_BMA <- function(x, dosevars) {
  if (is.null(x)) {
    return(NULL)
  }

  nm <- "included.realizations.BMA"
  check_integer(x, nm, minlen = 0, maxlen = 0, min = 1, max = length(dosevars))
  if (length(x) < 2) {
    stop(
      "ERROR: BMA requires at least two included exposure realizations. ",
      "With a single exposure realization, use RC."
    )
  }
  if (any(duplicated(x))) {
    stop(paste0("ERROR: ", nm, " contains duplicated values"))
  }

  NULL
}

check_M <- function(vars, data) {
  nm <- "M"
  check_vars(data, vars, nm, minlen = 0, maxlen = 0)
  for (v in vars) {
    vec <- data[, v, drop = TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2, binary = 1)
  }
  NULL
}

check_X <- function(vars, data) {
  nm <- "X"
  if (is.null(vars)) {
    return(NULL)
  }
  for (v in vars) {
    vec <- data[, v, drop = TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2)
  }
  NULL
}


warn_if_poorly_conditioned_X <- function(
  X_matrix,
  family,
  mean_sd_threshold = 100,
  sd_ratio_threshold = 1e4,
  kappa_threshold = 1e6
) {
  if (is.null(X_matrix) || !length(X_matrix)) {
    return(invisible(NULL))
  }

  X_matrix <- as.matrix(X_matrix)
  if (!nrow(X_matrix) || !ncol(X_matrix) || !all(is.finite(X_matrix))) {
    return(invisible(NULL))
  }

  sds <- apply(X_matrix, 2, sd)
  means <- colMeans(X_matrix)
  nonzero_sd <- is.finite(sds) & sds > 0
  details <- character()

  if (any(!nonzero_sd)) {
    details <- c(
      details,
      paste0(
        "near-zero variation in ",
        paste(colnames(X_matrix)[!nonzero_sd], collapse = ", ")
      )
    )
  }

  if (any(nonzero_sd)) {
    mean_sd_ratio <- abs(means[nonzero_sd]) / sds[nonzero_sd]
    large_mean <- mean_sd_ratio > mean_sd_threshold
    if (any(large_mean)) {
      details <- c(
        details,
        paste0(
          "large mean relative to standard deviation in ",
          paste(names(mean_sd_ratio)[large_mean], collapse = ", ")
        )
      )
    }

    sd_ratio <- max(sds[nonzero_sd]) / min(sds[nonzero_sd])
    if (is.finite(sd_ratio) && sd_ratio > sd_ratio_threshold) {
      details <- c(
        details,
        paste0(
          "standard deviations differ by a factor of ",
          signif(sd_ratio, 4)
        )
      )
    }
  }

  design <- if (
    family %in% c("gaussian", "binomial", "poisson", "multinomial")
  ) {
    cbind(`(Intercept)` = 1, X_matrix)
  } else {
    X_matrix
  }
  if (ncol(design) > 1 && nrow(design) >= ncol(design)) {
    design_kappa <- tryCatch(kappa(design), error = function(e) NA_real_)
    if (is.finite(design_kappa) && design_kappa > kappa_threshold) {
      details <- c(
        details,
        paste0("design condition number is ", signif(design_kappa, 4))
      )
    }
  }

  if (length(details)) {
    warning(
      paste0(
        "WARNING: right-hand-side covariates appear poorly scaled or ",
        "ill-conditioned (",
        paste(unique(details), collapse = "; "),
        "). This can make optimizer convergence unreliable. Consider ",
        "centering/scaling continuous covariates, such as calendar year, ",
        "before fitting."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


check_offset <- function(v, data) {
  if (!length(v)) {
    return(NULL)
  }
  nm <- "offset"
  check_vars(data, v, nm, minlen = 0, maxlen = 0)
  check_num_vec(data[, v, drop = TRUE], nm, nonneg = 1)

  NULL
}

check_setnr <- function(v, data) {
  nm <- "setnr"
  check_vars(data, v, nm, minlen = 1, maxlen = 1)
  check_num_vec(data[, v, drop = TRUE], nm, nonneg = 1, integer = 1)

  NULL
}

filter_clogit_sets <- function(data, status, setnr) {
  y <- data[, status, drop = TRUE]
  sets <- data[, setnr, drop = TRUE]
  set_sizes <- table(sets)
  cases_per_set <- tapply(y, sets, sum)

  multi_case_sets <- names(cases_per_set)[cases_per_set > 1]
  if (length(multi_case_sets)) {
    stop(paste0(
      "Conditional logistic regression currently requires exactly one case ",
      "per informative matched set; ",
      length(multi_case_sets),
      " matched set(s) contain more than one case."
    ))
  }

  size_one_sets <- names(set_sizes)[set_sizes == 1]
  no_case_sets <- names(cases_per_set)[cases_per_set == 0]
  drop_sets <- union(size_one_sets, no_case_sets)

  if (length(drop_sets)) {
    warning(paste0(
      "Excluding ",
      length(size_one_sets),
      " matched set(s) of size 1 and ",
      length(setdiff(no_case_sets, size_one_sets)),
      " additional matched set(s) with no cases from conditional logistic ",
      "regression."
    ))
    data <- data[!(sets %in% drop_sets), , drop = FALSE]
  }

  if (!nrow(data)) {
    stop(
      "No informative matched sets remain after excluding matched sets of ",
      "size 1 or with no cases."
    )
  }

  data
}

filter_prophaz_zero_followup <- function(data, entry, exit, status = NULL) {
  if (is.null(entry) || !length(entry)) {
    return(data)
  }

  entry_vec <- data[, entry, drop = TRUE]
  exit_vec <- data[, exit, drop = TRUE]
  zero_followup <- entry_vec == exit_vec
  zero_followup[is.na(zero_followup)] <- FALSE

  if (any(zero_followup)) {
    n_excluded <- sum(zero_followup)
    n_events <- if (!is.null(status) && length(status)) {
      sum(data[zero_followup, status, drop = TRUE] == 1, na.rm = TRUE)
    } else {
      NA_integer_
    }

    event_note <- if (is.na(n_events)) {
      ""
    } else {
      paste0(" Of these, ", n_events, " had status = 1.")
    }

    warning(paste0(
      "Excluding ",
      n_excluded,
      " row(s) with entry == exit since they have zero follow-up time.",
      event_note
    ))
    data <- data[!zero_followup, , drop = FALSE]
  }

  if (!nrow(data)) {
    stop(
      "No rows remain after excluding rows with entry == exit from ",
      "proportional hazards regression."
    )
  }

  data
}

check_entry_exit <- function(entry, exit, data) {
  nm1 <- "entry"
  nm2 <- "exit"
  check_vars(data, entry, nm1, minlen = 0, maxlen = 0)
  check_vars(data, exit, nm2, minlen = 1, maxlen = 1)
  vec2 <- data[, exit, drop = TRUE]
  check_num_vec(vec2, nm2)
  if (length(entry)) {
    vec1 <- data[, entry, drop = TRUE]
    check_num_vec(vec1, nm1)
    tmp <- vec1 > vec2
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) stop(paste0("ERROR: ", nm1, " > ", nm2, " for some values"))
  }

  NULL
}

check_methods <- function(x) {
  nm <- "methods"
  valid <- c("RC", "ERC", "MCML", "FMA", "BMA")
  ret <- check_char_vec(x, nm, valid = valid, def = "RC", len = 0)
  ret <- unique(ret)
  ret
}

check_deg <- function(x) {
  if (!length(x)) {
    return(2)
  }
  nm <- "deg"
  check_integer(x, nm, minlen = 1, maxlen = 1, min = 1, max = 2)
  x
}

check_inpar <- function(x, family, M, X, deg, multinom_levels = 0) {
  if (!length(x)) {
    return(NULL)
  }
  nm <- "inpar"
  if (family == "gaussian") {
    len <- 2 + length(X) + length(M) * deg + deg
  } else if (family %in% c("binomial", "poisson")) {
    len <- 1 + length(X) + length(M) * deg + deg
  } else if (family %in% c("prophaz", "clogit")) {
    len <- length(X) + length(M) * deg + deg
  } else if (family == "multinomial") {
    len <- (multinom_levels - 1) * (1 + length(X) + length(M) * deg + deg)
  } else {
    stop("ERROR")
  }
  check_num_vec(x, nm, binary = 0, nonneg = 0, integer = 0, len = len)
  x
}

check_future_chunk_size_FMA <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  nm <- "future.chunk.size.FMA"
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    stop(paste0("ERROR: ", nm, " must be a positive number or NULL"))
  }

  x
}


check_factor_vec <- function(x, nm, len = 0) {
  if (!is.factor(x)) {
    stop(paste0("ERROR: ", nm, " must be numeric"))
  }
  if (len && (len != length(x))) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector of length ", len))
  }

  if (length(levels(x)) < 3) {
    stop(paste0("ERROR: ", nm, " must have at least 3 levels"))
  }

  if (length(levels(x)) > length(unique(x))) {
    stop(paste0("ERROR: ", nm, " contains unused levels"))
  }

  NULL
}

check_num_vec <- function(x, nm, binary = 0, nonneg = 0, integer = 0, len = 0) {
  if (!is.numeric(x)) {
    stop(paste0("ERROR: ", nm, " must be numeric"))
  }
  if (len && (len != length(x))) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector of length ", len))
  }
  tmp <- !is.finite(x)
  if (any(tmp)) {
    stop(paste0("ERROR: ", nm, " must contain finite values"))
  }
  if (binary) {
    tmp <- !(x %in% 0:1)
    if (any(tmp)) {
      stop(paste0("ERROR: ", nm, " must contain binary (0 - 1) values"))
    }
  }
  if (nonneg) {
    tmp <- x < 0
    if (any(tmp)) {
      stop(paste0("ERROR: ", nm, " must contain non-negative values"))
    }
  }
  if (integer) {
    check_integer(x, nm, minlen = 0, maxlen = 0, min = NULL, max = NULL)
  }

  NULL
}

check_vars <- function(data, vars, nm, minlen = 0, maxlen = 0) {
  nv <- length(vars)
  if (minlen && (minlen == maxlen) && (nv != minlen)) {
    stop(paste0("ERROR: ", nm, " must have length ", minlen))
  }
  if (nv < minlen) {
    stop(paste0("ERROR: ", nm, " must have length >= ", minlen))
  }
  if (!nv) {
    return(NULL)
  }
  if (!is.vector(vars)) {
    stop(paste0(
      "ERROR: ",
      nm,
      " must be a vector of indices or variable names"
    ))
  }

  nc <- ncol(data)
  cx <- colnames(data)
  if (is.numeric(vars)) {
    check_integer(vars, nm, minlen = minlen, maxlen = nc, min = 1, max = nc)
  } else if (is.character(vars)) {
    check_char_vec(vars, nm, valid = cx, def = NULL)
  } else {
    stop(paste0(
      "ERROR: ",
      nm,
      " must be a vector of indices or variable names"
    ))
  }
  if (any(duplicated(vars))) {
    stop(paste0("ERROR: ", nm, " contains duplicated values"))
  }

  unique(vars)
}


check_char_vec <- function(x, nm, valid = NULL, def = NULL, len = 0) {
  n <- length(x)
  if (len && (n != len)) {
    stop(paste0("ERROR: ", nm, " must have length ", len))
  }
  if (!n) {
    return(def)
  }
  if (!is.character(x)) {
    stop(paste0("ERROR: ", nm, " must be character"))
  }
  if (length(valid)) {
    tmp <- !(x %in% valid)
    if (any(tmp)) {
      err <- getCharVecStr(x[tmp])
      stop(paste0("ERROR ", nm, " contains invalid values: ", err))
    }
  }
  x
}


required_vars <- function(m) {
  vars <- c(m$dosevars, model_X_vars(m), m$M_names)

  if (m$family %in% c("gaussian", "binomial", "poisson", "multinomial")) {
    vars <- c(vars, m$Y)
  }

  if (m$family == "poisson" && !is.null(m$offset)) {
    vars <- c(vars, m$offset)
  }

  if (m$family == "prophaz") {
    vars <- c(vars, m$status, m$exit)
    if (!is.null(m$entry)) vars <- c(vars, m$entry)
  }

  if (m$family == "clogit") {
    vars <- c(vars, m$status, m$setnr)
  }

  vars[!is.null(vars)]
}


model_X_vars <- function(m) {
  if (!is.null(m$X_names)) {
    return(m$X_names)
  }
  if (is.character(m$X)) {
    return(m$X)
  }
  character()
}


resolve_na_action <- function(na.action, env = parent.frame()) {
  if (is.null(na.action)) {
    na.action <- getOption("na.action")
  }
  if (is.null(na.action)) {
    na.action <- "na.omit"
  }

  if (is.character(na.action)) {
    if (length(na.action) != 1) {
      stop("ERROR: na.action must be a function or a single function name")
    }
    if (!exists(na.action, envir = env, mode = "function", inherits = TRUE)) {
      stop("ERROR: na.action function not found: ", na.action)
    }
    na.action <- get(na.action, envir = env, mode = "function", inherits = TRUE)
  }

  if (!is.function(na.action)) {
    stop("ERROR: na.action must be a function or a single function name")
  }

  na.action
}


model_na_vars <- function(
  family,
  dosevars,
  Y = NULL,
  status = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  entry = NULL,
  exit = NULL,
  setnr = NULL
) {
  vars <- c(dosevars, M, X)

  if (family %in% c("gaussian", "binomial", "poisson", "multinomial")) {
    vars <- c(vars, Y)
  }

  if (family == "poisson") {
    vars <- c(vars, offset)
  }

  if (family == "prophaz") {
    vars <- c(vars, status, exit, entry)
  }

  if (family == "clogit") {
    vars <- c(vars, status, setnr)
  }

  vars <- vars[!is.na(vars)]
  vars <- vars[nzchar(vars)]
  unique(vars)
}


model_na_vars_from_model <- function(m) {
  model_na_vars(
    family = m$family,
    dosevars = m$dosevars,
    Y = m$Y,
    status = m$status,
    M = m$M_names,
    X = model_X_vars(m),
    offset = m$offset,
    entry = m$entry,
    exit = m$exit,
    setnr = m$setnr
  )
}


apply_na_action_to_data <- function(data, vars, na.action) {
  vars <- unique(vars)
  vars <- vars[vars %in% colnames(data)]

  if (!length(vars)) {
    return(list(data = data, na.action = NULL))
  }

  na_frame <- data[, vars, drop = FALSE]
  row_id <- ".ameras_row_id"
  while (row_id %in% colnames(na_frame)) {
    row_id <- paste0(".", row_id)
  }
  na_frame[[row_id]] <- seq_len(nrow(data))

  acted <- na.action(na_frame)
  if (!is.data.frame(acted) || !row_id %in% colnames(acted)) {
    stop(
      "ERROR: na.action must return a data frame with rows preserved or removed"
    )
  }

  keep <- acted[[row_id]]
  if (!length(keep)) {
    stop("ERROR: no complete observations remain after applying na.action")
  }

  list(
    data = data[keep, , drop = FALSE],
    na.action = attr(acted, "na.action")
  )
}


check_integer <- function(
  x,
  nm,
  minlen = 1,
  maxlen = 0,
  min = NULL,
  max = NULL
) {
  n <- length(x)
  if (minlen && (minlen == maxlen) && (n != minlen)) {
    stop(paste0("ERROR: ", nm, " must have length ", minlen))
  }
  if (minlen && (n < minlen)) {
    stop(paste0("ERROR: ", nm, " must have length >= ", minlen))
  }
  if (!is.numeric(x)) {
    stop(paste0("ERROR: ", nm, " must be integer"))
  }
  if (any(!is.finite(x))) {
    stop(paste0("ERROR: ", nm, " must be integer"))
  }
  if (any(x != floor(x))) {
    stop(paste0("ERROR: ", nm, " must be integer"))
  }
  if (length(min) && any(x < min)) {
    stop(paste0("ERROR: ", nm, " must be >= ", min))
  }
  if (length(max) && any(x > max)) {
    stop(paste0("ERROR: ", nm, " must be <= ", max))
  }

  NULL
}


check_string <- function(obj, valid, parm) {
  # obj:   A character string (length 1)
  # valid: Character vector of valid values
  # parm:  The name of the argument being checked

  errFlag <- 0

  # Check for errors
  if (!isString(obj)) {
    errFlag <- 1
  }
  if (!errFlag) {
    obj <- trimws(obj)
    if (!(obj %in% valid)) errFlag <- 1
  }
  if (errFlag) {
    msg <- getCharVecStr(valid)
    msg <- paste0("ERROR: ", parm, " contains the invalid values ", msg)
    stop(msg)
  }

  obj
} # END: check.string

getVarNumbers <- function(vars, data) {
  if (!length(vars)) {
    return(NULL)
  }
  if (is.numeric(vars)) {
    return(vars)
  }
  cx <- colnames(data)
  ret <- match(vars, cx)
  ret
}


getCharVecStr <- function(x, sep = ",") {
  ret <- paste0("'", x, "'")
  ret <- paste0(ret, collapse = sep)
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
