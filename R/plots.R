ecdfplot <- function(
  data,
  dosevars,
  xlab = "Dose",
  ylab = "Cumulative distribution",
  show.mean = TRUE,
  log.xaxis = TRUE
) {
  check_pkgs(c("dplyr", "ggplot2", "tidyr", "scales", "patchwork"))

  dose <- curve_id <- row_id <- row_mean <- NULL # To suppress R CMD CHECK notes about undefined global variables

  dosemat <- data[, dosevars, drop = FALSE]

  p1 <- dosemat |>
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
    ggplot2::ggtitle("Curves corresponding to realizations")
  if (show.mean) {
    means1 <- dosemat |>
      dplyr::mutate(
        row_mean = rowMeans(dplyr::across(dplyr::everything()), na.rm = TRUE)
      )
    p1 <- p1 +
      ggplot2::stat_ecdf(
        data = means1,
        ggplot2::aes(x = row_mean),
        color = "black",
        linewidth = 1.2
      )
  }
  if (log.xaxis) {
    p1 <- p1 + ggplot2::scale_x_log10(labels = scales::comma)
  }

  p2 <- as.data.frame(t(dosemat)) |>
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
    ggplot2::ggtitle("Curves corresponding to individuals")
  if (show.mean) {
    means2 <- as.data.frame(t(dosemat)) |>
      dplyr::mutate(
        row_mean = rowMeans(dplyr::across(dplyr::everything()), na.rm = TRUE)
      )
    p2 <- p2 +
      ggplot2::stat_ecdf(
        data = means2,
        ggplot2::aes(x = row_mean),
        color = "black",
        linewidth = 1.2
      )
  }

  if (log.xaxis) {
    p2 <- p2 + ggplot2::scale_x_log10(labels = scales::comma)
  }
  patchwork::wrap_plots(p1, p2)
}
