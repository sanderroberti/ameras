



ecdfplot <- function(data, dosevars, xlab="Dose", ylab="Cumulative distribution"){

  check_pkgs(c("dplyr","ggplot2", "tidyr", "scales"))

  dose <- curve_id <- row_id <- NULL # To suppress R CMD CHECK notes about undefined global variables

  dosemat <- data[, dosevars, drop=FALSE ]

  dosemat |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -row_id,
      names_to = "curve_id",
      values_to = "dose"
    ) |>
    ggplot2::ggplot() +
    ggplot2::stat_ecdf(
      ggplot2::aes(x = dose, group = curve_id),
      color = "grey",
      alpha = 0.2
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_x_log10(labels = scales::comma)
}

