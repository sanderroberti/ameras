make_method_parnames <- function(
  family,
  data,
  Y = NULL,
  X = NULL,
  M = NULL,
  deg = 1,
  doseRRmod = NULL,
  include_sigma = identical(family, "gaussian"),
  prophaz_numints = NULL,
  modifier_info = NULL
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
  parnames <- c(
    parnames,
    x_names,
    modifier_dose_names(
      doseRRmod = doseRRmod,
      deg = deg,
      data = data,
      M = M,
      modifier_info = modifier_info
    )
  )

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
