test_that("ecdfplot returns a patchwork object with and without mean curves", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("scales")

  data("data", package = "ameras")
  plot_data <- data[1:8, ]
  dosevars <- paste0("V", 1:3)

  p_with_means <- ecdfplot(
    plot_data,
    dosevars,
    xlab = "Custom dose",
    ylab = "Custom cumulative",
    show.mean = TRUE,
    log.xaxis = FALSE
  )

  p_without_means <- ecdfplot(
    plot_data,
    dosevars,
    show.mean = FALSE,
    log.xaxis = TRUE
  )

  expect_s3_class(p_with_means, "patchwork")
  expect_s3_class(p_without_means, "patchwork")
})
