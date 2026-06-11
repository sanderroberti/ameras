make_method_parnames <- function(
  family,
  data,
  Y = NULL,
  X = NULL,
  M = NULL,
  deg = 1,
  doseRRmod = NULL,
  include_sigma = identical(family, "gaussian"),
  prophaz_numints = NULL
) {
  # Most families have an intercept-like parameter; conditional logistic and
  # proportional hazards methods do not. Keep that distinction explicit because
  # it differs from the rest of the dose-response naming rules.
  parnames <- if (family %in% c("prophaz", "clogit")) {
    NULL
  } else {
    "(Intercept)"
  }

  x_names <- names(data[, X, drop = FALSE])
  m_names <- names(data[, M, drop = FALSE])

  if (identical(doseRRmod, "LINEXP")) {
    parnames <- c(
      parnames,
      x_names,
      c("dose_linear", "dose_exponential")
    )
    if (!is.null(M)) {
      parnames <- c(parnames, paste0("dose_linear:", m_names))
      parnames <- c(parnames, paste0("dose_exponential:", m_names))
    }
  } else {
    parnames <- c(
      parnames,
      x_names,
      c("dose", "dose_squared")[1:deg]
    )
    if (!is.null(M)) {
      parnames <- c(parnames, paste0("dose:", m_names))
      if (deg == 2) {
        parnames <- c(parnames, paste0("dose_squared:", m_names))
      }
    }
  }

  if (include_sigma) {
    parnames <- c(parnames, "sigma")
  }
  if (!is.null(prophaz_numints)) {
    parnames <- c(parnames, paste0("h0[", seq_len(prophaz_numints), "]"))
  }

  if (identical(family, "multinomial")) {
    # Multinomial fits repeat the base names for each modeled response level.
    # The final factor level is the reference level and keeps no parameter block.
    mylv <- levels(data[, Y])
    mylv <- mylv[-length(mylv)]
    parnames <- do.call(
      "c",
      lapply(mylv, function(lv) paste0("(", lv, ")_", parnames))
    )
  }

  parnames
}
