test_that("method timing combines fit and CI phases", {
  fit <- ameras:::new_runtime_phase(cpu = 1.2, elapsed = 5)
  ci <- ameras:::new_runtime_phase(cpu = 0.3, elapsed = 2)

  timing <- ameras:::new_method_timing(fit = fit, ci = ci)

  expect_equal(timing$fit$cpu, 1.2)
  expect_equal(timing$ci$cpu, 0.3)
  expect_equal(timing$total$cpu, 1.5)
  expect_equal(timing$total$elapsed, 7)
  expect_identical(ameras:::format_runtime(timing$total$cpu), "1.5 seconds")
})

test_that("CI timing replaces previous CI timing instead of accumulating", {
  method_fit <- list(
    runtime = "1 seconds",
    timing = ameras:::new_method_timing(
      fit = ameras:::new_runtime_phase(cpu = 1, elapsed = 10),
      ci = ameras:::new_runtime_phase(cpu = 4, elapsed = 40)
    )
  )

  updated <- ameras:::set_ci_timing(
    method_fit,
    ameras:::new_runtime_phase(cpu = 0.25, elapsed = 3)
  )

  expect_equal(updated$timing$fit$cpu, 1)
  expect_equal(updated$timing$ci$cpu, 0.25)
  expect_equal(updated$timing$total$cpu, 1.25)
  expect_identical(updated$runtime, "1.2 seconds")
})

test_that("summary uses structured timing when available", {
  object <- ameras:::new_amerasfit(list(
    call = quote(ameras(Y ~ dose(D), data = dat)),
    num.rows = 1,
    num.realizations = 1,
    CI.computed = FALSE,
    RC = list(
      coefficients = c(a = 1),
      sd = c(a = 0.1),
      vcov = matrix(0.01, dimnames = list("a", "a")),
      runtime = "99 seconds",
      timing = ameras:::new_method_timing(
        fit = ameras:::new_runtime_phase(cpu = 1.2, elapsed = 10),
        ci = ameras:::new_runtime_phase(cpu = 0.3, elapsed = 2)
      )
    )
  ))

  summ <- summary(object)

  expect_named(summ$runtime_table, c("Method", "Fit", "CI", "Total"))
  expect_equal(summ$runtime_table$Fit, 1.2)
  expect_equal(summ$runtime_table$CI, 0.3)
  expect_equal(summ$runtime_table$Total, 1.5)
  expect_equal(summ$total_runtime_seconds, 1.5)
})

test_that("summary falls back to legacy runtime strings", {
  object <- ameras:::new_amerasfit(list(
    call = quote(ameras(Y ~ dose(D), data = dat)),
    num.rows = 1,
    num.realizations = 1,
    CI.computed = FALSE,
    RC = list(
      coefficients = c(a = 1),
      sd = c(a = 0.1),
      vcov = matrix(0.01, dimnames = list("a", "a")),
      runtime = "2 seconds"
    )
  ))

  summ <- summary(object)

  expect_named(summ$runtime_table, c("Method", "Runtime"))
  expect_equal(summ$runtime_table$Runtime, 2)
  expect_equal(summ$total_runtime_seconds, 2)
})

test_that("confint records CI timing for sample-based intervals", {
  object <- ameras:::new_amerasfit(list(
    call = quote(ameras(Y ~ dose(D), data = dat, methods = "FMA")),
    num.rows = 3,
    num.realizations = 1,
    CI.computed = FALSE,
    FMA = list(
      coefficients = c(a = 1),
      sd = c(a = 0.1),
      vcov = matrix(0.01, dimnames = list("a", "a")),
      samples = data.frame(a = c(0.8, 1, 1.2)),
      runtime = "1 seconds",
      timing = ameras:::new_method_timing(
        fit = ameras:::new_runtime_phase(cpu = 1, elapsed = 2)
      )
    )
  ))

  out <- suppressMessages(confint(object, type = "percentile", print = FALSE))

  expect_true(out$CI.computed)
  expect_equal(out$FMA$timing$fit$cpu, 1)
  expect_true(out$FMA$timing$ci$cpu >= 0)
  expect_equal(
    out$FMA$timing$total$cpu,
    out$FMA$timing$fit$cpu + out$FMA$timing$ci$cpu
  )
  expect_identical(
    out$FMA$runtime,
    ameras:::format_runtime(out$FMA$timing$total$cpu)
  )
})
