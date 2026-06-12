new_runtime_phase <- function(cpu = 0, elapsed = 0) {
  list(
    cpu = as.numeric(cpu),
    elapsed = as.numeric(elapsed)
  )
}

start_runtime_timer <- function() {
  proc.time()
}

stop_runtime_timer <- function(start) {
  elapsed <- proc.time() - start
  cpu_names <- c("user.self", "sys.self", "user.child", "sys.child")
  cpu <- sum(
    as.numeric(elapsed[intersect(cpu_names, names(elapsed))]),
    na.rm = TRUE
  )

  new_runtime_phase(
    cpu = cpu,
    elapsed = as.numeric(elapsed["elapsed"])
  )
}

add_runtime_phases <- function(...) {
  phases <- list(...)

  new_runtime_phase(
    cpu = sum(vapply(phases, function(x) x$cpu %||% 0, numeric(1))),
    elapsed = sum(vapply(phases, function(x) x$elapsed %||% 0, numeric(1)))
  )
}

new_method_timing <- function(
  fit = new_runtime_phase(),
  ci = new_runtime_phase()
) {
  list(
    fit = fit,
    ci = ci,
    total = add_runtime_phases(fit, ci)
  )
}

format_runtime <- function(seconds) {
  paste(round(as.numeric(seconds), 1), "seconds")
}

runtime_seconds_from_string <- function(runtime) {
  if (is.null(runtime) || is.na(runtime)) {
    return(NA_real_)
  }

  as.numeric(strsplit(runtime, " seconds")[[1]])
}

timing_from_legacy_runtime <- function(runtime) {
  seconds <- runtime_seconds_from_string(runtime)

  new_method_timing(
    fit = new_runtime_phase(cpu = seconds, elapsed = seconds)
  )
}

set_ci_timing <- function(method_fit, ci_timing) {
  timing <- method_fit$timing %||% timing_from_legacy_runtime(method_fit$runtime)
  method_fit$timing <- new_method_timing(
    fit = timing$fit,
    ci = ci_timing
  )
  method_fit$runtime <- format_runtime(method_fit$timing$total$cpu)
  method_fit
}

method_runtime_row <- function(method_fit, method) {
  if (!is.null(method_fit$timing)) {
    return(data.frame(
      Method = method,
      Fit = method_fit$timing$fit$cpu,
      CI = method_fit$timing$ci$cpu,
      Total = method_fit$timing$total$cpu,
      row.names = NULL
    ))
  }

  data.frame(
    Method = method,
    Runtime = runtime_seconds_from_string(method_fit$runtime),
    row.names = NULL
  )
}

runtime_table_from_methods <- function(object0) {
  rows <- lapply(seq_along(object0), function(i) {
    method_runtime_row(object0[[i]], names(object0)[i])
  })

  use_structured <- any(vapply(
    rows,
    function(x) "Total" %in% names(x),
    logical(1)
  ))

  if (use_structured) {
    rows <- lapply(rows, function(x) {
      if ("Total" %in% names(x)) {
        return(x)
      }

      data.frame(
        Method = x$Method,
        Fit = x$Runtime,
        CI = 0,
        Total = x$Runtime,
        row.names = NULL
      )
    })
  }

  do.call("rbind", rows)
}

total_runtime_seconds <- function(runtime_table) {
  if ("Total" %in% names(runtime_table)) {
    return(sum(runtime_table$Total))
  }

  sum(runtime_table$Runtime)
}
