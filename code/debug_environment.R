#!/usr/bin/env Rscript
# Diagnose the existing report's map cache and first environmental request.
# This is a standalone diagnostic, not a replacement for pipeline code.

debug_log <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%dT%H:%M:%S%z]"), " [environment-debug] ",
      paste0(..., collapse = ""), "\n", sep = "", file = stderr())
}

debug_probe <- function(state, name, action) {
  debug_log("START ", name)
  started <- proc.time()[["elapsed"]]
  failure <- NULL
  value <- tryCatch(withCallingHandlers(action(), warning = function(w) {
    debug_log("WARNING ", name, ": ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { failure <<- e; NULL })
  elapsed <- proc.time()[["elapsed"]] - started
  ok <- is.null(failure)
  detail <- if (ok) "" else conditionMessage(failure)
  debug_log(if (ok) "PASS " else "FAIL ", name, " elapsed=",
            round(elapsed, 2), "s", if (!ok) paste0(" | ", detail) else "")
  if (!ok && !is.null(conditionCall(failure))) {
    debug_log("ERROR CALL: ", paste(deparse(conditionCall(failure)), collapse = " "))
  }
  state$results[[length(state$results) + 1L]] <- data.frame(
    probe = name, status = if (ok) "PASS" else "FAIL", seconds = elapsed,
    detail = gsub("[\t\r\n]", " ", detail), stringsAsFactors = FALSE)
  if (!is.null(state$out)) {
    write.table(do.call(rbind, state$results), file.path(state$out, "results.tsv"),
                sep = "\t", row.names = FALSE, quote = TRUE)
  }
  list(ok = ok, value = value, error = if (ok) NULL else failure)
}

read_rmd_functions <- function(path, names) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  fence <- strrep(intToUtf8(96), 3)
  inside <- FALSE
  code <- character()
  for (line in lines) {
    if (grepl(paste0("^", fence, "\\{r([[:space:]]|\\})"), line)) {
      inside <- TRUE
    } else if (trimws(line) == fence) {
      inside <- FALSE
    } else if (inside) code <- c(code, line)
  }
  parsed <- parse(text = code)
  found <- setNames(vector("list", length(names)), names)
  visit <- function(x) {
    if (is.call(x) && identical(x[[1]], as.name("<-")) &&
        is.symbol(x[[2]]) && as.character(x[[2]]) %in% names &&
        is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) {
      key <- as.character(x[[2]])
      if (!is.null(found[[key]])) stop("Multiple function definitions for ", key)
      found[[key]] <<- x[[3]]
    }
    if (is.call(x) || is.expression(x) || is.pairlist(x)) {
      for (j in seq_along(x)) {
        if (!identical(x[[j]], quote(expr = ))) visit(x[[j]])
      }
    }
  }
  visit(parsed)
  missing <- names[vapply(found, is.null, logical(1))]
  if (length(missing)) stop("Rmd functions not found: ", paste(missing, collapse = ", "))
  found
}

rewrite_nc_create <- function(expr) {
  count <- 0L
  walk <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(x[[1]], quote(ncdf4::nc_create)) ||
        identical(x[[1]], as.name("nc_create"))) {
      x[[1]] <- as.name(".debug_nc_create")
      count <<- count + 1L
    }
    for (j in seq_along(x)) {
      if (!identical(x[[j]], quote(expr = ))) x[[j]] <- walk(x[[j]])
    }
    x
  }
  result <- walk(expr)
  list(expr = result, count = count)
}

clone_ncdf_writer <- function(writer, force_v4 = NULL,
                              create = function(...) ncdf4::nc_create(...)) {
  replacement <- rewrite_nc_create(body(writer))
  if (replacement$count != 1L) {
    stop("Expected one nc_create call in the installed writer, found ", replacement$count,
         ". Original writer remains unchanged.")
  }
  original_environment <- environment(writer)
  shim <- new.env(parent = original_environment)
  shim$.debug_nc_create <- function(...) {
    args <- list(...)
    args$verbose <- TRUE
    if (!is.null(force_v4)) args$force_v4 <- force_v4
    debug_log("nc_create: force_v4=",
              if (is.null(args$force_v4)) "package default" else args$force_v4)
    do.call(create, args)
  }
  body(writer) <- replacement$expr
  environment(writer) <- shim
  writer
}

select_debug_coords <- function(coords, sample_order, roi = rep(NA_real_, 4)) {
  required <- c("Sample ID", "Latitude", "Longitude")
  if (!all(required %in% names(coords))) stop("Coordinate columns are not standardized")
  coords$Latitude <- suppressWarnings(as.numeric(coords$Latitude))
  coords$Longitude <- suppressWarnings(as.numeric(coords$Longitude))
  matched <- match(sample_order, coords$`Sample ID`)
  debug_log("VCF samples=", length(sample_order), "; unmatched coordinate IDs=", sum(is.na(matched)))
  coords <- coords[matched, , drop = FALSE]
  coords <- coords[is.finite(coords$Latitude) & is.finite(coords$Longitude), , drop = FALSE]
  if (length(roi) != 4L) stop("ROI must contain lat_min, lat_max, lon_min, lon_max")
  if (all(is.finite(roi))) {
    if (roi[1] > roi[2] || roi[3] > roi[4]) stop("ROI minima exceed maxima")
    coords <- coords[coords$Latitude >= roi[1] & coords$Latitude <= roi[2] &
                       coords$Longitude >= roi[3] & coords$Longitude <= roi[4], , drop = FALSE]
  }
  if (!nrow(coords)) stop("No matched samples with usable coordinates remain")
  if (any(abs(coords$Latitude) > 90 | abs(coords$Longitude) > 180)) {
    stop("Coordinates outside WGS84 latitude/longitude bounds")
  }
  coords
}

debug_region <- function(coords) {
  longitude <- range(coords$Longitude)
  latitude <- range(coords$Latitude)
  px <- diff(longitude) * 0.05
  py <- diff(latitude) * 0.05
  unname(c(longitude[1] - px, latitude[1] - py, longitude[2] + px, latitude[2] + py))
}

small_debug_region <- function(region) {
  if (length(region) != 4L || any(!is.finite(region)) ||
      region[1] >= region[3] || region[2] >= region[4]) stop("Invalid request region")
  mid <- c(mean(region[c(1, 3)]), mean(region[c(2, 4)]))
  half <- pmin(c(region[3] - region[1], region[4] - region[2]) / 2, 0.5)
  unname(c(mid - half, mid + half))
}

debug_storage <- function(path) {
  debug_log("Filesystem: ", normalizePath(path, mustWork = TRUE),
            "; writable=", file.access(path, 2) == 0)
  if (nzchar(Sys.which("df"))) {
    for (option in c("-Pk", "-Pi")) {
      output <- system2(Sys.which("df"), c(option, shQuote(path)), stdout = TRUE, stderr = TRUE)
      cat(paste(output, collapse = "\n"), "\n")
    }
  }
}

debug_tiny_netcdf <- function(directory, force_v4) {
  file <- tempfile("tiny-ncdf-", tmpdir = directory, fileext = ".nc")
  on.exit(unlink(file), add = TRUE)
  dimension <- ncdf4::ncdim_def("x", "1", 1:3)
  variable <- ncdf4::ncvar_def("value", "1", list(dimension), -999, prec = "double")
  nc <- ncdf4::nc_create(file, list(variable), force_v4 = force_v4, verbose = TRUE)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)
  ncdf4::ncvar_put(nc, "value", c(1, 2, 3))
  ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(file)
  stopifnot(identical(as.numeric(ncdf4::ncvar_get(nc, "value")), c(1, 2, 3)))
  debug_log("Tiny NetCDF round trip: ", nc$format)
  invisible(TRUE)
}

debug_tiny_vectors <- function(directory) {
  paths <- c(tempfile("tiny-sf-", tmpdir = directory, fileext = ".gpkg"),
             tempfile("tiny-terra-", tmpdir = directory, fileext = ".gpkg"))
  on.exit(unlink(paths), add = TRUE)
  point <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326))
  sf::st_write(point, paths[1], quiet = FALSE)
  stopifnot(nrow(sf::st_read(paths[1], quiet = FALSE)) == 1L)
  terra::writeVector(terra::vect(point), paths[2], overwrite = TRUE)
  stopifnot(nrow(terra::vect(paths[2])) == 1L)
  invisible(TRUE)
}

debug_raw_summary <- function(raw) {
  debug_log("Returned class: ", paste(class(raw), collapse = ", "),
            "; object.size MiB=", round(as.numeric(object.size(raw)) / 1024^2, 2))
  debug_log("Dimensions: ", paste(names(dim(raw)), dim(raw), sep = "=", collapse = ", "))
  print(stars::st_dimensions(raw))
  for (name in names(raw)) {
    array <- raw[[name]]
    bytes <- as.double(length(array)) * 8
    debug_log("Variable ", name, ": type=", typeof(array), "; dimensions=",
              paste(dim(array), collapse = "x"), "; values=", length(array),
              "; estimated double output GiB=", signif(bytes / 1024^3, 5),
              "; above roughly 2GiB=", bytes >= 2^31)
  }
}

debug_arguments <- function(args) {
  result <- list(config = "config/config.yaml", rmd = "code/2_Herring_Wingen.Rmd",
                 out = file.path("debug", paste0("environment-", format(Sys.time(), "%Y%m%d-%H%M%S"),
                                                 "-", Sys.getpid())),
                 gene_vcf = NULL, small_only = FALSE, skip_maps = FALSE)
  j <- 1L
  while (j <= length(args)) {
    key <- args[j]
    if (key %in% c("--small-only", "--skip-maps")) {
      result[[gsub("-", "_", substring(key, 3))]] <- TRUE
      j <- j + 1L
    } else if (key %in% c("--config", "--rmd", "--out", "--gene-vcf")) {
      if (j == length(args)) stop("Missing value for ", key)
      result[[gsub("-", "_", substring(key, 3))]] <- args[j + 1L]
      j <- j + 2L
    } else stop("Unknown option: ", key)
  }
  result
}

debug_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- debug_arguments(args)
  dir.create(options$out, recursive = TRUE, showWarnings = FALSE)
  options$out <- normalizePath(options$out, mustWork = TRUE)
  state <- new.env(parent = emptyenv())
  state$results <- list()
  state$out <- options$out
  debug_log("START standalone diagnostic; PID=", Sys.getpid(), "; output=", options$out)
  debug_log("R tempdir=", tempdir(), "; TMPDIR=", Sys.getenv("TMPDIR"))
  packages <- c("yaml", "sf", "terra", "ncdf4", "stars", "CopernicusMarine", "rnaturalearth")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Required packages missing from this runtime: ", paste(absent, collapse = ", "),
                           ". Use the same container as the failing reports.")
  for (package in packages) debug_log(package, "=", as.character(utils::packageVersion(package)))
  debug_log("NetCDF library: ", ncdf4::nc_version())
  print(sf::sf_extSoftVersion())
  writeLines(capture.output(sessionInfo()), file.path(options$out, "session-info.txt"))

  config <- yaml::read_yaml(options$config)
  root <- Sys.getenv("ROOT", unset = getwd())
  resolve <- function(path) if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
  params <- config$wingen_params
  setting <- function(name, default) if (is.null(params[[name]])) default else params[[name]]
  preset <- setting("ENV_PRESET", "global")
  presets_path <- Sys.getenv("ENV_PRESETS_FILE", unset = file.path(root, "config", "env_presets.yaml"))
  presets <- yaml::read_yaml(presets_path)$presets
  chosen <- presets[[preset]]
  if (is.null(chosen)) stop("Preset not found: ", preset)
  years <- as.numeric(c(setting("ENV_BASELINE_START_YEAR", 2019), setting("ENV_BASELINE_END_YEAR", 2022)))
  if (!is.null(chosen$physics$first_year)) years[1] <- max(years[1], chosen$physics$first_year)
  if (!is.null(chosen$physics$last_year)) years[2] <- min(years[2], chosen$physics$last_year)
  if (any(!is.finite(years)) || years[1] > years[2]) stop("Invalid physics baseline")
  depth <- as.numeric(setting("ENV_DEPTH_M", 0))
  if (!is.finite(depth) || depth < 0) stop("Invalid physics depth")
  results <- resolve(config$results_dir)
  ne_cache <- file.path(results, "temp", "ne_cache")
  env_cache <- file.path(results, "temp", "env_cache")
  debug_log("Config=", normalizePath(options$config), "; Rmd=", normalizePath(options$rmd))
  versions <- tools::md5sum(c(options$config, options$rmd, presets_path))
  write.table(data.frame(path = names(versions), md5 = unname(versions)),
              file.path(options$out, "input-files.tsv"), sep = "\t", row.names = FALSE)
  debug_log("Preset=", preset, "; physics product=", chosen$physics$product,
            "; layer=", chosen$physics$layer, "; salinity variable=", chosen$physics$salinity_var,
            "; baseline=", paste(years, collapse = "-"), "; depth=", depth)
  for (path in c(ne_cache, env_cache)) {
    debug_log("Existing cache entries in ", path, ":")
    info <- file.info(list.files(path, full.names = TRUE))
    if (nrow(info)) print(info[, c("size", "mtime"), drop = FALSE])
  }

  # These functions are parsed from the user's Rmd, never from a copied implementation.
  source_functions <- read_rmd_functions(options$rmd, c("standardize_coord_columns", "download_env_series"))
  source_environment <- new.env(parent = globalenv())
  source_environment$input_coordinate_file_path <- resolve(config$sample_coordinates_file)
  standardize <- eval(source_functions$standardize_coord_columns, source_environment)
  gene_vcf <- options$gene_vcf
  if (is.null(gene_vcf)) {
    files <- sort(list.files(file.path(results, "temp", "filtered_vcf"),
                             pattern = "_filtered\\.vcf(\\.gz)?$", full.names = TRUE, recursive = TRUE))
    if (!length(files)) stop("No completed filtered gene VCF found; supply --gene-vcf PATH")
    gene_vcf <- files[1]
  }
  gene_vcf <- resolve(gene_vcf)
  if (!file.exists(gene_vcf)) stop("Filtered gene VCF not found: ", gene_vcf)
  if (!nzchar(Sys.which("bcftools"))) stop("bcftools is required only to read the existing gene's sample names")
  samples <- system2(Sys.which("bcftools"), c("query", "-l", shQuote(gene_vcf)), stdout = TRUE)
  if ((!is.null(attr(samples, "status")) && attr(samples, "status") != 0L) || !length(samples)) {
    stop("Could not read sample names from the existing gene VCF")
  }
  debug_log("Using completed gene VCF: ", gene_vcf, "; no variants are regenerated")
  coordinate_path <- source_environment$input_coordinate_file_path
  extension <- tolower(tools::file_ext(coordinate_path))
  coordinates <- if (extension %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl required for coordinate workbook")
    readxl::read_excel(coordinate_path)
  } else if (extension == "csv") {
    if (!requireNamespace("readr", quietly = TRUE)) stop("readr required to match the report's CSV reader")
    readr::read_delim(coordinate_path, show_col_types = FALSE)
  } else stop("Unsupported coordinate file extension: ", extension)
  roi <- vapply(c("lat_min", "lat_max", "lon_min", "lon_max"), function(key) {
    if (is.null(config[[key]]) || !nzchar(as.character(config[[key]]))) NA_real_ else as.numeric(config[[key]])
  }, numeric(1))
  coordinates <- select_debug_coords(standardize(coordinates), samples, roi)
  region <- debug_region(coordinates)
  debug_log("Retained coordinate rows=", nrow(coordinates),
            "; padded request xmin,ymin,xmax,ymax=", paste(region, collapse = ", "))
  request <- list(preset = preset, gene_vcf = gene_vcf, product = chosen$physics$product,
                  layer = chosen$physics$layer, variable = chosen$physics$salinity_var,
                  years = years, depth = depth, region = region, roi = roi)
  dput(request, file.path(options$out, "request.R"))

  # The unique probe subdirectories never replace or delete existing cache files.
  directories <- list(R_tempdir = tempdir(), diagnostic_output = options$out,
                      Natural_Earth_cache = ne_cache, environmental_cache = env_cache)
  for (label in names(directories)) {
    parent <- directories[[label]]
    debug_probe(state, paste0("filesystem_", label), function() {
      dir.create(parent, recursive = TRUE, showWarnings = FALSE)
      debug_storage(parent)
      probe_dir <- tempfile(".environment-debug-", tmpdir = parent)
      if (!dir.create(probe_dir)) stop("Cannot create unique probe directory in ", parent)
      on.exit(unlink(probe_dir, recursive = TRUE), add = TRUE)
      for (force in c(FALSE, TRUE)) {
        debug_probe(state, paste0("tiny_netcdf_", label, "_", if (force) "v4" else "classic"),
                    function() debug_tiny_netcdf(probe_dir, force))
      }
      debug_probe(state, paste0("tiny_geopackage_", label), function() debug_tiny_vectors(probe_dir))
      TRUE
    })
  }

  if (!options$skip_maps) {
    for (type in c("land", "ocean")) {
      returnclass <- if (type == "land") "sf" else "sv"
      debug_probe(state, paste0("shared_map_read_", type), function() {
        object <- rnaturalearth::ne_load(type = type, category = "physical", scale = 10,
                                       destdir = ne_cache, returnclass = returnclass)
        debug_log("Shared map rows=", nrow(object))
        if (!nrow(object)) stop("Shared map has zero features")
        TRUE
      })
      debug_probe(state, paste0("private_map_download_", type), function() {
        private <- tempfile(paste0("private-", type, "-"), tmpdir = options$out)
        if (!dir.create(private)) stop("Cannot create private map directory")
        # Retain the private download for inspection; shared cache is not modified.
        debug_log("Fresh private map download: ", private)
        rnaturalearth::ne_download(type = type, category = "physical", scale = 10,
                                  destdir = private, load = FALSE, returnclass = returnclass)
        object <- rnaturalearth::ne_load(type = type, category = "physical", scale = 10,
                                       destdir = private, returnclass = returnclass)
        debug_log("Private map rows=", nrow(object))
        if (!nrow(object)) stop("Fresh private map has zero features")
        TRUE
      })
    }
  }

  writer <- getExportedValue("CopernicusMarine", "cms_write_ncdf")
  writeLines(deparse(writer), file.path(options$out, "installed-netcdf-writer.R"))
  # Only the clone is modified. Installed package namespaces remain untouched.
  writer_v4 <- clone_ncdf_writer(writer, force_v4 = TRUE)
  probe_environment <- new.env(parent = globalenv())
  eval_call <- eval(source_functions$download_env_series, probe_environment)
  run_request <- function(label, bounds, timerange) {
    raw <- NULL
    paths <- character()
    stage <- "before download"
    writer_failed <- FALSE
    on.exit({ unlink(paths); rm(raw); gc(verbose = FALSE) }, add = TRUE)
    probe_environment$env_region <- bounds
    probe_environment$env_verticalrange <- c(-depth, -depth)
    probe_environment$cms_download_subset <- function(...) {
      stage <<- "download subset"
      request_args <- list(...)
      request_args$timerange <- timerange
      debug_log("Request ", label, ": product=", chosen$physics$product,
                "; region=", paste(bounds, collapse = ","), "; time=", paste(timerange, collapse = "/"))
      raw <<- do.call(CopernicusMarine::cms_download_subset, request_args)
      debug_log("Download returned for ", label)
      stage <<- "inspect returned array"
      debug_raw_summary(raw)
      raw
    }
    probe_environment$cms_write_ncdf <- function(x, file, ...) {
      stage <<- "original NetCDF write"
      paths <<- c(paths, file)
      debug_log("Original package writer destination=", file)
      debug_storage(dirname(file))
      # Standalone Rscript exposes the package's Rprintf output, unlike failed knit chunks.
      value <- tryCatch(CopernicusMarine::cms_write_ncdf(x, file, ...), error = function(e) {
        writer_failed <<- TRUE
        stop(e)
      })
      debug_log("Original NetCDF write completed; bytes=", file.info(file)$size)
      stage <<- "open original NetCDF with terra::rast"
      invisible(value)
    }
    result <- debug_probe(state, paste0(label, "_original_writer"), function() {
      raster <- eval_call(chosen$physics$product, chosen$physics$layer,
                          chosen$physics$salinity_var, years, c(-depth, -depth))
      debug_log("Original raster: rows=", terra::nrow(raster), "; columns=", terra::ncol(raster),
                "; layers=", terra::nlyr(raster))
      TRUE
    })
    if (!result$ok) debug_log("Failed stage for ", label, ": ", stage)
    if (writer_failed && !is.null(raw)) {
      debug_probe(state, paste0(label, "_same_data_netcdf4"), function() {
        file <- tempfile(paste0(label, "-v4-"), fileext = ".nc")
        paths <<- c(paths, file)
        debug_log("Same object, same R temp directory, diagnostic NetCDF4 destination=", file)
        writer_v4(raw, file)
        debug_log("NetCDF4 write completed; bytes=", file.info(file)$size,
                  "; opening with terra::rast")
        raster <- terra::rast(file)
        debug_log("NetCDF4 raster: rows=", terra::nrow(raster), "; columns=", terra::ncol(raster),
                  "; layers=", terra::nlyr(raster))
        TRUE
      })
    }
    invisible(NULL)
  }
  run_request("small_request", small_debug_region(region),
              c(sprintf("%d-01-01", as.integer(years[1])), sprintf("%d-01-31", as.integer(years[1]))))
  if (!options$small_only) {
    run_request("full_request", region,
                c(sprintf("%d-01-01", as.integer(years[1])), sprintf("%d-12-31", as.integer(years[2]))))
  }
  summary <- do.call(rbind, state$results)
  print(summary, row.names = FALSE)
  failures <- sum(summary$status == "FAIL")
  debug_log("DONE probes=", nrow(summary), "; failed probes=", failures,
            "; see results.tsv and debug.log. A failed probe is diagnostic evidence.")
  if (failures) 1L else 0L
}

if (sys.nframe() == 0L) {
  status <- tryCatch(debug_main(), error = function(e) {
    debug_log("DIAGNOSTIC SETUP ERROR: ", conditionMessage(e))
    2L
  })
  quit(save = "no", status = status)
}
