grid_timings <- new.env(parent = emptyenv())
grid_timings$records <- data.frame(stage = character(), engine = character(),
  elapsed_seconds = numeric(), status = character(), stringsAsFactors = FALSE)
grid_timed_call <- function(stage, engine, fn, ...) {
  started <- proc.time()[["elapsed"]]
  completed <- FALSE
  on.exit({
    elapsed <- proc.time()[["elapsed"]] - started
    status <- if (completed) "ok" else "error"
    grid_timings$records <- rbind(grid_timings$records,
      data.frame(stage = stage, engine = engine, elapsed_seconds = elapsed,
                 status = status, stringsAsFactors = FALSE))
    grid_log(sprintf("[grid-timing] stage=%s; engine=%s; elapsed=%.3fs; status=%s",
                     stage, engine, elapsed, status))
  }, add = TRUE)
  result <- withVisible(fn(...))
  completed <- TRUE
  if (result$visible) result$value else invisible(result$value)
}
grid_timed_engine <- function(fn, stage, engine) {
  force(fn); force(stage); force(engine)
  function(...) grid_timed_call(stage, engine, fn, ...)
}
grid_timing_sum <- function(stage, engine) {
  rows <- grid_timings$records$stage == stage &
          grid_timings$records$engine == engine & grid_timings$records$status == "ok"
  if (!any(rows)) return(NA_real_)
  sum(grid_timings$records$elapsed_seconds[rows])
}

grid_assert_raster_equal <- function(expected, actual, label) {
  if (!terra::compareGeom(expected, actual, stopOnError = FALSE) ||
      !identical(names(expected), names(actual))) {
    stop(label, ": raster geometry or layer names differ")
  }
  # Read rows in bounded batches; avoid two full-grid validation matrices.
  nr <- terra::nrow(expected)
  rows <- max(1, floor(50000 / terra::ncol(expected)))
  terra::readStart(expected)
  terra::readStart(actual)
  on.exit({ terra::readStop(expected); terra::readStop(actual) }, add = TRUE)
  for (row in seq.int(1, nr, by = rows)) {
    nrows <- min(rows, nr - row + 1)
    a <- terra::readValues(expected, row = row, nrows = nrows, mat = TRUE)
    b <- terra::readValues(actual, row = row, nrows = nrows, mat = TRUE)
    if (!identical(is.na(a), is.na(b)) ||
        !identical(a[!is.na(a)], b[!is.na(b)])) {
      stop(label, ": values or missing-cell masks differ at row ", row)
    }
  }
  grid_log("[grid-validate] PASS ", label,
          ": exact nonmissing values, missing masks, geometry and layer names")
  invisible(TRUE)
}

grid_checked_window <- function(optimized, native) {
  force(optimized); force(native)
  function(...) {
    before <- get(".Random.seed", envir = .GlobalEnv)
    actual <- grid_timed_call("windows", "optimized", optimized, ...)
    after <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", after, envir = .GlobalEnv), add = TRUE)
    assign(".Random.seed", before, envir = .GlobalEnv)
    expected <- grid_timed_call("windows", "native_reference", native, ...)
    if (!identical(after, get(".Random.seed", envir = .GlobalEnv)))
      stop("Window validation: post-calculation RNG state differs")
    grid_assert_raster_equal(expected, actual, "moving windows")
    actual
  }
}

grid_checked_kriging <- function(optimized, domain) {
  force(optimized); force(domain)
  function(...) {
    args <- list(...)
    want_model <- isTRUE(args$model_output)
    args$model_output <- TRUE
    actual <- grid_timed_call("kriging", "optimized", function() do.call(optimized, args))
    expected <- grid_timed_call("kriging", "native_reference", function() do.call(wingen::wkrig_gd, args))
    if (!identical(expected$model, actual$model) ||
        !identical(expected$variogram, actual$variogram))
      stop("Kriging validation: variogram or fitted model differs")
    expected$raster <- terra::mask(expected$raster, domain, touches = TRUE)
    grid_assert_raster_equal(expected$raster, actual$raster, "kriging after domain mask")
    if (want_model) actual else actual$raster
  }
}


# Diagnostics profile calls by default. Native comparisons need an explicit
# GRID_VALIDATE=TRUE as well and can take hours on a large rectangular grid.
decorate_wingen_functions <- function(functions, domain, optimized) {
  value <- tolower(Sys.getenv("GRID_VALIDATE", unset = "FALSE"))
  if (!value %in% c("true", "false", "1", "0")) stop("GRID_VALIDATE must be TRUE or FALSE")
  validate <- value %in% c("true", "1") && optimized
  grid_log("[diagnostics] native comparison=", validate)
  if (validate) {
    functions$window_gd <- grid_checked_window(functions$window_gd, wingen::window_gd)
    functions$window_general <- grid_checked_window(functions$window_general, wingen::window_general)
    functions$wkrig_gd <- grid_checked_kriging(functions$wkrig_gd, domain)
  } else {
    engine <- if (optimized) "optimized" else "native"
    functions$window_gd <- grid_timed_engine(functions$window_gd, "windows", engine)
    functions$window_general <- grid_timed_engine(functions$window_general, "windows", engine)
    functions$wkrig_gd <- grid_timed_engine(functions$wkrig_gd, "kriging", engine)
  }
  attr(functions, "validation") <- validate
  functions
}

report_sampling_geometry <- function(sites, indices, distances, spacing, threshold_km) {
  nearest <- apply(distances, 1L, function(d) min(c(d[d > 0], Inf)))
  grid_log("[sampling] individuals=", length(indices), "; exact sites=", nrow(sites),
           "; spacing groups at ", threshold_km, " km=", length(unique(spacing$groups)))
  grid_log("[sampling] median individual-weighted nearest positive spacing=",
           round(stats::median(nearest[indices]) / 1e3, 3), " km; clustered spacing=",
           round(spacing$reach_km, 3), " km")
}

write_diagnostic_timings <- function(directory, prefix) {
  utils::write.table(grid_timings$records,
    file.path(directory, paste0(prefix, "grid_engine_timings.tsv")),
    sep = "\t", row.names = FALSE, quote = TRUE)
}

# Optional console run: Rscript code/help_functions/wingen_diagnostics.R code/2_Herring_Wingen.Rmd [logfile]
run_diagnostics <- function(args) {
  if (length(args) < 1L || length(args) > 2L)
    stop("Supply the Rmd path and optionally a console-log filename")
  rmd <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
  log <- if (length(args) == 2L) args[2L] else
    file.path(dirname(rmd), paste0(format(Sys.Date(), "%Y-%m-%d"), "_wingen_diagnostics_console.txt"))
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  child <- tempfile(fileext = ".R")
  on.exit(unlink(child), add = TRUE)
  # Tee knitr's captured text as well as process output, retaining report output.
  writeLines(c(
    'Sys.setenv(WINGEN_DIAGNOSTICS = "TRUE", R_EVALUATE_BYPASS_MESSAGES = "true")',
    'knitr::knit_hooks$set(evaluate = local({',
    '  evaluate <- knitr::knit_hooks$get("evaluate")',
    '  function(..., output_handler) {',
    '    text <- output_handler$text',
    '    output_handler$text <- function(x) { cat(x, file = stderr()); text(x) }',
    '    evaluate(..., output_handler = output_handler)',
    '  }',
    '}))',
    paste0("rmarkdown::render(", encodeString(rmd, quote = '"'), ")")), child)
  exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(exe, shQuote(child), stdout = log, stderr = log)
  cat("Diagnostic console log: ", normalizePath(log, winslash = "/"), "\n", sep = "")
  status
}
if (!interactive() && sys.nframe() == 0L)
  quit(save = "no", status = run_diagnostics(commandArgs(trailingOnly = TRUE)))
