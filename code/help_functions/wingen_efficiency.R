grid_log <- function(...) cat(..., "\n", sep = "")

# Use one working-grid cell count for BOTH stages, so a threshold comparison
# cannot accidentally leave windows native while kriging is optimized.
grid_use_optimized <- function(cells, mode = "auto", min_cells = 5000) {
  if (length(cells) != 1L || !is.numeric(cells) || !is.finite(cells) ||
      cells < 1 || cells != floor(cells))
    stop("Working grid cell count must be a positive whole number")
  if (length(mode) != 1L || is.na(mode))
    stop("GRID_EFFICIENCY must be auto, TRUE or FALSE")
  mode <- tolower(as.character(mode))
  if (!mode %in% c("auto", "true", "false", "1", "0"))
    stop("GRID_EFFICIENCY must be auto, TRUE or FALSE")
  if (length(min_cells) != 1L || !is.numeric(min_cells) || is.na(min_cells) ||
      min_cells < 0 || (is.finite(min_cells) && min_cells != floor(min_cells)))
    stop("GRID_MIN_CELLS must be a nonnegative whole number or Inf")
  enabled <- if (mode == "auto") cells >= min_cells else mode %in% c("true", "1")
  grid_log("[grid-efficiency] working_grid_cells=", cells,
          "; threshold=", min_cells, "; mode=", mode,
          "; window+kriging=", if (enabled) "optimized" else "native")
  enabled
}


# Rectangular Wingen acceleration for the SAME calculation raster.
# No cell size, extent, origin, mask or sample coordinate is changed. Only
# neighborhoods known to contain enough samples are evaluated. This is tied to
# audited Wingen 2.2.0 / furrr 0.4.0 internals and fails closed on other versions.

window_wingen_internal <- function(name) {
  value <- get0(name, envir = asNamespace("wingen"), inherits = FALSE)
  if (!is.function(value)) stop("Installed Wingen lacks helper '", name, "'")
  value
}

sparse_window_membership <- function(lyr, coords, wdim) {
  if (!inherits(lyr, "SpatRaster") || terra::nlyr(lyr) != 1L)
    stop("lyr must be a single-layer SpatRaster")
  if (isTRUE(terra::is.lonlat(lyr)))
    stop("Sparse rectangular membership requires a projected raster")
  if (!is.numeric(wdim) || !length(wdim) %in% c(1L, 2L) ||
      anyNA(wdim) || any(!is.finite(wdim)) ||
      any(wdim < 1 | wdim != floor(wdim) | wdim %% 2 != 1))
    stop("wdim must contain one or two positive odd whole numbers")
  if (length(wdim) == 1L) wdim <- rep(wdim, 2L)
  wdim <- as.integer(wdim)
  if (all(wdim == 1L)) stop("wdim must span more than one cell")

  coord_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]
  nr <- terra::nrow(lyr)
  nc <- terra::ncol(lyr)
  mask_na <- is.na(terra::values(lyr, mat = FALSE))
  count_values <- rep(0L, terra::ncell(lyr))
  count_values[mask_na] <- NA_integer_
  valid_samples <- which(!is.na(coord_cells))
  if (length(valid_samples)) {
    sample_rc <- terra::rowColFromCell(lyr, coord_cells[valid_samples])
    offsets <- expand.grid(
      row = seq.int(-(wdim[1L] %/% 2L), wdim[1L] %/% 2L),
      col = seq.int(-(wdim[2L] %/% 2L), wdim[2L] %/% 2L))
    # Reverse lookup: a sample contributes to neighboring unmasked focal cells.
    # A sample ON a masked cell still contributes to nearby valid focal cells.
    n_offsets <- nrow(offsets)
    center_rows <- rep(sample_rc[, 1L], each = n_offsets) +
      rep(offsets$row, times = length(valid_samples))
    center_cols <- rep(sample_rc[, 2L], each = n_offsets) +
      rep(offsets$col, times = length(valid_samples))
    individuals <- rep(valid_samples, each = n_offsets)
    in_grid <- center_rows >= 1L & center_rows <= nr &
      center_cols >= 1L & center_cols <= nc
    center_ids <- (center_rows[in_grid] - 1L) * nc + center_cols[in_grid]
    individuals <- individuals[in_grid]
    unmasked <- !mask_na[center_ids]
    center_ids <- center_ids[unmasked]
    individuals <- individuals[unmasked]
    # Native get_adj() returns which(...), i.e. original individual order.
    ord <- order(center_ids, individuals)
    cells <- sort(unique(center_ids))
    members <- unname(split(individuals[ord],
                            factor(center_ids[ord], levels = cells)))
    count_values[cells] <- lengths(members)
  } else {
    cells <- integer()
    members <- list()
  }
  sample_count <- terra::setValues(terra::rast(lyr), count_values)
  names(sample_count) <- "sample_count"
  list(cells = cells, members = members, sample_count = sample_count,
       coord_cells = coord_cells, wdim = wdim)
}

# Reproduce furrr::make_seeds(TRUE, n) at ORIGINAL raster cell IDs. Retain only
# selected seeds, avoiding a list of seven-integer seeds for every grid cell.
# Stream advancement still costs O(largest selected cell ID); it is necessary
# to keep the same rarefaction draws as native Wingen with the same set.seed().
window_sparse_seeds <- function(n, cells) {
  stopifnot(length(n) == 1L, is.finite(n), n >= 0, n == floor(n),
            is.numeric(cells), all(is.finite(cells)),
            all(cells == floor(cells)), all(cells >= 1 & cells <= n),
            !anyDuplicated(cells))
  ns <- asNamespace("furrr")
  rng_start <- get("as_lecyer_cmrg_seed_from_base_r_seed", ns)
  rng_next <- get("next_random_seed", ns)
  rng_set <- get("set_random_seed", ns)
  seed <- rng_start()
  after_make_seeds <- rng_next()
  on.exit(rng_set(after_make_seeds), add = TRUE)
  out <- vector("list", length(cells))
  if (!length(cells)) return(out)
  ord <- order(cells)
  sorted <- cells[ord]
  j <- 1L
  for (i in seq_len(max(sorted))) {
    if (i == sorted[[j]]) {
      out[[ord[[j]]]] <- parallel::nextRNGSubStream(seed)
      j <- j + 1L
      if (j > length(sorted)) break
    }
    seed <- parallel::nextRNGStream(seed)
  }
  out
}

window_sparse_map <- function(n, cells, fn, ...) {
  # furrr's normal finalization must run even for an empty selected set.
  if (!length(cells))
    return(furrr::future_map(cells, fn, ...,
      .options = furrr::furrr_options(seed = TRUE, packages = c("wingen", "terra"))))
  seeds <- window_sparse_seeds(n, cells)
  furrr::future_map(cells, fn, ...,
    .options = furrr::furrr_options(seed = seeds, packages = c("wingen", "terra")))
}

make_grid_efficient_windows <- function(prepared = NULL) {
  invisible(lapply(c("wingen", "furrr", "terra"), loadNamespace))
  if (as.character(utils::packageVersion("wingen")) != "2.2.0" ||
      as.character(utils::packageVersion("furrr")) != "0.4.0")
    stop("Grid-efficient windows were audited for Wingen 2.2.0 / furrr 0.4.0; ",
         "use native Wingen and revalidate before changing these version guards")
  needed <- c("window_gd", "window_gd_stats", "window_general", "run_general",
              "vcf_check", "check_data", "convert_vcf", "return_stat",
              "rarify_helper", "sample_gd", "rm_duplicate_sample_count")
  invisible(lapply(needed, window_wingen_internal))
  dispatch <- new.env(parent = asNamespace("wingen"))
  cache <- new.env(parent = emptyenv())
  cache$geometry <- prepared
  native_run <- window_wingen_internal("run_general")
  dispatch$run_general <- function(x, lyr, coords,
                                  coord_cells = NULL, nmat = NULL,
                                  distmat = NULL, maxdist = NULL,
                                  stat, rarify, rarify_n, rarify_nit, min_n,
                                  fun, L, rarify_alleles, sig, ...) {
    # Preserve native behavior for argument combinations outside this planar,
    # positive-minimum rectangular fast path (e.g. geographic wraparound).
    if (is.null(nmat) || !is.null(distmat) || isTRUE(terra::is.lonlat(lyr)) ||
        length(min_n) != 1L || !is.finite(min_n) || min_n < 1)
      return(native_run(x, lyr, coords, coord_cells, nmat, distmat, maxdist,
                        stat, rarify, rarify_n, rarify_nit, min_n,
                        fun, L, rarify_alleles, sig, ...))
    started <- proc.time()[[3L]]
    corrected <- window_wingen_internal("check_data")(x, coords = coords)
    x <- corrected$x
    coords <- corrected$coords
    current_cells <- terra::extract(lyr, coords, cell = TRUE)[, "cell"]
    geometry <- cache$geometry
    can_reuse <- !is.null(geometry) &&
      identical(as.numeric(current_cells), as.numeric(geometry$coord_cells)) &&
      identical(as.integer(dim(nmat)), geometry$wdim) &&
      terra::compareGeom(lyr, geometry$sample_count, stopOnError = FALSE) &&
      identical(is.na(terra::values(lyr, mat = FALSE)),
                is.na(terra::values(geometry$sample_count, mat = FALSE)))
    if (!can_reuse) {
      geometry <- sparse_window_membership(lyr, coords, dim(nmat))
      cache$geometry <- geometry
    }
    counts <- terra::values(geometry$sample_count, mat = FALSE)
    if (!any(counts >= min_n, na.rm = TRUE)) {
      # Native preview_count() tests masked centers too. Preserve its unusual
      # error behavior while avoiding its all-cell loop for ordinary runs.
      check_grid <- terra::setValues(terra::rast(lyr), 1)
      check_counts <- sparse_window_membership(check_grid, coords, dim(nmat))$sample_count
      if (!any(terra::values(check_counts, mat = FALSE) >= min_n, na.rm = TRUE))
        stop("Minimum sample size (min_n) is not met for any window across this raster")
    }
    if (is.character(stat) && !is.null(L) && stat == "pi" && L == "nvariants")
      L <- ncol(x)
    stat_function <- window_wingen_internal("return_stat")(stat = stat, ...)
    if (is.null(rarify_n)) rarify_n <- min_n
    required <- if (rarify) max(min_n, rarify_n) else min_n
    selected <- which(lengths(geometry$members) >= required)
    cells <- geometry$cells[selected]
    members <- geometry$members[selected]
    names(members) <- as.character(cells)
    grid_log(sprintf(paste0("[window-gd] START grid=%s cells; unmasked=%s; ",
                           "sample-near=%s; eligible=%s; required_individuals=%s; ",
                           "geometry_reused=%s; workers=%s"),
      format(terra::ncell(lyr), big.mark = ",", scientific = FALSE),
      format(sum(!is.na(counts)), big.mark = ",", scientific = FALSE),
      length(geometry$cells), length(cells), required, can_reuse, future::nbrOfWorkers()))
    results <- window_sparse_map(terra::ncell(lyr), cells, function(i) {
      sub <- members[[as.character(i)]]
      if (rarify) {
        gd <- window_wingen_internal("rarify_helper")(
          x, sub, rarify_n, rarify_nit, stat_function, fun,
          L = L, rarify_alleles = rarify_alleles, sig = sig)
      } else {
        gd <- window_wingen_internal("sample_gd")(
          x, sub, stat_function, L = L, rarify_alleles = rarify_alleles, sig = sig)
      }
      if (is.null(names(gd))) names(gd) <- "custom"
      if (all(is.na(gd))) numeric() else gd
    })
    stat_names <- unique(unlist(lapply(results, names), use.names = FALSE))
    layers <- lapply(stat_names, function(name) {
      values <- rep(NA_real_, terra::ncell(lyr))
      values[cells] <- vapply(results, function(gd) {
        if (name %in% names(gd)) as.numeric(gd[[name]]) else NA_real_
      }, numeric(1))
      layer <- terra::setValues(terra::rast(lyr), values)
      names(layer) <- name
      layer
    })
    result <- terra::rast(c(layers, list(geometry$sample_count)))
    if (!length(stat_names))
      warning("All values NA, only a sample count layer will be returned")
    grid_log(sprintf("[window-gd] DONE evaluated=%s cells; elapsed=%.3fs",
                    length(cells), proc.time()[[3L]] - started))
    result
  }
  # Keep native public validation, aggregation, conversion, edge cropping and
  # stat-vector/layer ordering by cloning into a private dispatch environment.
  for (name in c("window_gd", "window_gd_stats", "window_general")) {
    entry <- window_wingen_internal(name)
    environment(entry) <- dispatch
    assign(name, entry, envir = dispatch)
  }
  window_gd <- function(gen, coords, lyr, stat = "pi", wdim = 3, fact = 0,
                        rarify = FALSE, rarify_n = NULL, rarify_nit = 5, min_n = 2,
                        fun = mean, L = "nvariants", rarify_alleles = TRUE, sig = 0.05,
                        crop_edges = FALSE, ...) {
    # Native window_gd_stats repeats these scans for each statistic. Run them
    # once per public call; conversion and statistics remain native. The call
    # environment is transient and never cached between genes.
    vcf <- window_wingen_internal("vcf_check")(gen)
    corrected <- window_wingen_internal("check_data")(vcf, coords)
    call_dispatch <- new.env(parent = dispatch)
    call_dispatch$vcf_check <- function(x) corrected$vcf
    call_dispatch$check_data <- function(x, coords = NULL, distmat = NULL) {
      if (inherits(x, "vcfR")) return(corrected)
      window_wingen_internal("check_data")(x, coords, distmat)
    }
    for (name in c("window_gd", "window_gd_stats", "window_general")) {
      entry <- window_wingen_internal(name)
      environment(entry) <- call_dispatch
      assign(name, entry, envir = call_dispatch)
    }
    call_dispatch$window_gd(gen, coords, lyr, stat, wdim, fact, rarify,
                           rarify_n, rarify_nit, min_n, fun, L,
                           rarify_alleles, sig, crop_edges, ...)
  }
  list(window_gd = window_gd, window_general = dispatch$window_general)
}

# Compatibility name used by the earlier reviewed helper.
make_sparse_wingen <- make_grid_efficient_windows


# Bound prediction memory without changing Wingen's variogram fit or training
# observations. The installed package is never modified.
grid_predict_blocks <- function(model, grd, prediction_crs, domain = NULL,
                                chunk_cells = 50000L, filename = "") {
  if (length(chunk_cells) != 1L || !is.finite(chunk_cells) ||
      chunk_cells < 1 || chunk_cells != floor(chunk_cells)) {
    stop("KRIGING_CHUNK_CELLS must be a positive whole number")
  }
  template <- terra::rast(grd)[[1L]]
  if (!is.null(domain) && !inherits(domain, "SpatVector")) {
    stop("Prediction domain must be a SpatVector polygon")
  }
  # Recompute at FINAL resolution: disaggregating a coarse NA mask would lose
  # some coastal cells that the original final touches=TRUE mask retains.
  eligible <- NULL
  if (!is.null(domain)) {
    eligible <- terra::mask(terra::init(terra::rast(template), 1), domain,
      touches = TRUE, filename = tempfile("wingen-domain-", fileext = ".tif"),
      overwrite = TRUE, datatype = "INT1U")
    terra::readStart(eligible)
  }
  if (!nzchar(filename)) filename <- tempfile("wingen-prediction-", fileext = ".tif")
  out <- terra::rast(template)
  # FLT8S keeps R's double precision when streaming to disk.
  terra::writeStart(out, filename, overwrite = TRUE,
                    wopt = list(datatype = "FLT8S"))
  writing <- TRUE
  on.exit({
    if (writing) try(terra::writeStop(out), silent = TRUE)
    if (!is.null(eligible)) {
      try(terra::readStop(eligible), silent = TRUE)
      path <- terra::sources(eligible)
      if (length(path)) unlink(path)
    }
  }, add = TRUE)
  nr <- terra::nrow(template)
  nc <- terra::ncol(template)
  rows_per_block <- max(1, floor(chunk_cells / nc))
  predicted <- 0
  for (row in seq.int(1, nr, by = rows_per_block)) {
    nrows <- min(rows_per_block, nr - row + 1)
    nvalues <- nrows * nc
    vals <- rep(NA_real_, nvalues)
    keep <- if (is.null(eligible)) seq_len(nvalues) else
      which(!is.na(terra::readValues(eligible, row = row, nrows = nrows,
                                   mat = FALSE)))
    # Even exceptionally wide rasters never create >chunk_cells sf points.
    if (length(keep)) {
      for (start in seq.int(1, length(keep), by = chunk_cells)) {
        ix <- keep[seq.int(start, min(length(keep), start + chunk_cells - 1))]
        cells <- (row - 1) * nc + ix
        xy <- terra::xyFromCell(template, cells)
        points <- sf::st_as_sf(data.frame(x = xy[, 1L], y = xy[, 2L]),
                              coords = c("x", "y"), crs = prediction_crs)
        predictions <- terra::predict(model, newdata = points, debug.level = 0)
        if (nrow(predictions) != length(ix) ||
            !"var1.pred" %in% names(predictions)) {
          stop("Unexpected prediction output; cannot preserve cell alignment")
        }
        vals[ix] <- predictions[["var1.pred"]]
        predicted <- predicted + length(ix)
      }
    }
    terra::writeValues(out, vals, start = row, nrows = nrows)
  }
  out <- terra::writeStop(out)
  writing <- FALSE
  grid_log(sprintf("[kriging-grid] predicted=%s of %s cells; sf batch limit=%s; FLT8S disk raster",
                  predicted, terra::ncell(template), chunk_cells))
  out
}

make_grid_efficient_kriging <- function(domain = NULL, chunk_cells = 50000L) {
  if (as.character(utils::packageVersion("wingen")) != "2.2.0") {
    stop("Grid-efficient kriging was validated with Wingen 2.2.0; set GRID_EFFICIENCY=FALSE for another version")
  }
  # Keep all native input checks, weighting, model selection, parameters and
  # model_output handling. Replace ONLY the four full-grid prediction steps.
  original <- getExportedValue("wingen", "wkrig_gd")
  exprs <- as.list(body(original))
  targets <- c("grd_df", "grd_sf", "krig_pred", "r_pred")
  positions <- vapply(targets, function(name) {
    hits <- which(vapply(exprs, function(e) is.call(e) &&
      identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name(name)), logical(1)))
    if (length(hits) != 1L) stop("Wingen prediction implementation changed: ", name)
    hits
  }, integer(1))
  if (!identical(unname(diff(positions)), rep(1L, 3L))) {
    stop("Wingen prediction steps changed; refusing an unverified rewrite")
  }
  expected_calls <- c("terra::as.data.frame", "sf::st_as_sf",
                      "terra::predict", "terra::rasterize")
  for (j in seq_along(positions)) {
    rhs <- exprs[[positions[j]]][[3L]]
    if (!is.call(rhs) || !identical(paste(deparse(rhs[[1L]]), collapse = ""),
                                   expected_calls[j])) {
      stop("Wingen prediction implementation changed: ", targets[j])
    }
  }
  replacement <- quote(r_pred <- grid_predict_blocks(krig_model, grd,
    terra::crs(r), domain = .grid_domain, chunk_cells = .grid_chunk_cells))
  exprs <- append(exprs[-positions], list(replacement), after = positions[1L] - 1L)
  optimized <- original
  body(optimized) <- as.call(exprs)
  dispatch <- new.env(parent = environment(original))
  dispatch$grid_predict_blocks <- grid_predict_blocks
  dispatch$.grid_domain <- domain
  dispatch$.grid_chunk_cells <- chunk_cells
  environment(optimized) <- dispatch
  optimized
}


# Configure local Wingen-compatible bindings; never modify its namespace.
configure_wingen <- function(lyr, domain, min_cells = 5000, chunk_cells = 50000, diagnostics = FALSE) {
  if (length(chunk_cells) != 1L || !is.finite(chunk_cells) ||
      chunk_cells < 1 || chunk_cells != floor(chunk_cells))
    stop("KRIGING_CHUNK_CELLS must be a positive whole number")
  enabled <- grid_use_optimized(terra::ncell(lyr), "auto", min_cells)
  compatible <- as.character(utils::packageVersion("wingen")) == "2.2.0" &&
                as.character(utils::packageVersion("furrr")) == "0.4.0"
  if (enabled && !compatible) {
    grid_log("[grid-efficiency] package versions differ from the tested adapters; using native Wingen")
    enabled <- FALSE
  }
  functions <- list(window_gd = wingen::window_gd,
                    window_general = wingen::window_general,
                    wkrig_gd = wingen::wkrig_gd)
  if (enabled) {
    windows <- make_grid_efficient_windows()
    functions$window_gd <- windows$window_gd
    functions$window_general <- windows$window_general
    functions$wkrig_gd <- make_grid_efficient_kriging(domain, chunk_cells)
  }
  decorator <- get0("decorate_wingen_functions", envir = environment(configure_wingen), inherits = FALSE)
  if (diagnostics) {
    if (!is.function(decorator)) stop("Reload the load-wingen-helpers chunk with diagnostic_run=TRUE")
    functions <- decorator(functions, domain, enabled)
  }
  attr(functions, "optimized") <- enabled
  functions
}

# Exact numeric matching: no coordinate rounding or genotype deduplication.
exact_site_index <- function(individuals, sites) {
  index <- integer(nrow(individuals))
  for (i in seq_len(nrow(sites)))
    index[individuals$Latitude == sites$Latitude[i] &
          individuals$Longitude == sites$Longitude[i]] <- i
  if (any(index == 0L)) stop("An individual has no exact sampling coordinate")
  index
}

# Only needed for the alternative spacing-based interpolation cutoff or a
# diagnostic report. Location-statistic/plot groups use their own shared IDs.
sample_spacing <- function(sites, tree, threshold_km) {
  groups <- if (nrow(sites) > 1L) stats::cutree(tree, h = threshold_km) else 1L
  centroids <- as.data.frame(t(vapply(sort(unique(groups)), function(i)
    c(lon = mean(sites$Longitude[groups == i]), lat = mean(sites$Latitude[groups == i])),
    numeric(2))))
  points <- terra::vect(centroids, geom = c("lon", "lat"), crs = "EPSG:4326")
  distances <- as.matrix(terra::distance(points))
  nearest <- apply(distances, 1L, function(d) min(c(d[d > 0], Inf)))
  list(reach_km = stats::median(nearest) / 1e3, groups = groups,
       centroids = centroids, nearest_m = nearest)
}
