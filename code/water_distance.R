# water_distance.R
# Shortest route distance (km) from a START - one or several points and/or boxed regions - to many sample locations, travelling only
# through water ("ocean" domain) or only over land ("land" domain). Never as the crow flies.
#
# Method: rasterise the coastline polygons onto a SQUARE-cell metric grid, build a graph of the passable cells
# (16 neighbours; an edge exists only if the straight segment between the two cell centres stays inside passable
# cells), and run ONE Dijkstra search per start. That single search gives the distance to every cell, so the cost is
# independent of the number of sample locations (seconds for the whole Baltic at 1 km).
# A start POINT is the passable cell nearest to it. A start BOX is a lon/lat rectangle - e.g. the HELCOM line that marks the entrance
# to the Baltic Sea - and stands for every passable cell it touches: they all lie at distance 0 (one virtual vertex joined to each
# of them by a zero-length edge), so the distance of a site is the shortest route to ANY cell of the box.
#
# Why not terra::gridDist() directly? Same idea, but (tested on terra 1.7-65; the current source has the same rule)
#   * a diagonal step only needs the target cell to be non-NA -> it slips between two land cells that touch at a corner,
#   * it only has 8 neighbours: open-water distances come out up to +8 % too long (mean ~+5 %; 16 neighbours: +2.7 % / +1.4 %),
#   * it gives distances only, no routes to inspect.
# Do NOT run it on the Mollweide grid used for wingen: Mollweide is equal-AREA, straight-line distances in it are
# off by -16 % .. +14 % across the Baltic (EPSG:3035 / a local LAEA: < 0.5 %).
# The mask rule matters as much as the algorithm: with Natural Earth 10 m at 1 km the Hel spit only stays a barrier if
# "any land touching a cell makes it land" (min_cover = 0.99); with a cell-centre rule it leaks (see test_water_distance.R).
#
# Needs: terra, igraph (+ jsonlite for the cache helpers). Everything is called with `pkg::` and nothing is attached,
# so sourcing this file inside the wingen Rmd cannot mask any function the Rmd already uses.

WATER_DIST_VERSION <- "5"     # part of the cache key: bump it whenever the algorithm, its defaults or the table columns change

.wd_need <- function(pkgs = c("terra", "igraph")) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
    stop("water_distance.R needs the R package '", p, "' (add it to the container / conda environment).", call. = FALSE)
}

## ---- 1. graph on a logical passable-cell matrix ------------------------------------------------------------
## W: logical matrix, TRUE = passable, row 1 = northern-most row. res: cell size in metres.
## corner_safe: an edge exists only if EVERY cell the centre-to-centre segment crosses is passable. For a diagonal step the two cells
##   flanking it must both be passable (strict, deliberately: relaxing it to 'not both blocked' shortened real Baltic routes by only 0.1-0.3 %).
grid_graph <- function(W, res, connectivity = 16, corner_safe = TRUE) {
  .wd_need("igraph")
  stopifnot(connectivity %in% c(8, 16))
  nr <- nrow(W); nc <- ncol(W)
  Wl <- W; Wl[is.na(Wl)] <- FALSE
  vid <- matrix(0L, nr, nc); vid[Wl] <- seq_len(sum(Wl))        # vertex i  <->  which(Wl)[i] (column-major)
  offs <- list(c(0, 1), c(1, 0), c(1, 1), c(1, -1))
  if (connectivity == 16) offs <- c(offs, list(c(1, 2), c(2, 1), c(1, -2), c(2, -1)))   # knight moves
  fr <- to <- wt <- vector("list", length(offs))
  for (k in seq_along(offs)) {
    dr <- offs[[k]][1]; dc <- offs[[k]][2]
    ri <- 1:(nr - dr); ci <- if (dc >= 0) 1:(nc - dc) else (1 - dc):nc
    ok <- Wl[ri, ci, drop = FALSE] & Wl[ri + dr, ci + dc, drop = FALSE]
    if (corner_safe) {
      if (dr == 1 && dc ==  1) ok <- ok & Wl[ri + 1, ci, drop = FALSE] & Wl[ri, ci + 1, drop = FALSE]
      if (dr == 1 && dc == -1) ok <- ok & Wl[ri + 1, ci, drop = FALSE] & Wl[ri, ci - 1, drop = FALSE]
      if (dr == 1 && dc ==  2) ok <- ok & Wl[ri, ci + 1, drop = FALSE] & Wl[ri + 1, ci + 1, drop = FALSE]
      if (dr == 2 && dc ==  1) ok <- ok & Wl[ri + 1, ci, drop = FALSE] & Wl[ri + 1, ci + 1, drop = FALSE]
      if (dr == 1 && dc == -2) ok <- ok & Wl[ri, ci - 1, drop = FALSE] & Wl[ri + 1, ci - 1, drop = FALSE]
      if (dr == 2 && dc == -1) ok <- ok & Wl[ri + 1, ci, drop = FALSE] & Wl[ri + 1, ci - 1, drop = FALSE]
    }
    fr[[k]] <- vid[ri, ci, drop = FALSE][ok]
    to[[k]] <- vid[ri + dr, ci + dc, drop = FALSE][ok]
    wt[[k]] <- rep(res * sqrt(dr^2 + dc^2), length(fr[[k]]))
  }
  g <- igraph::make_graph(as.vector(rbind(unlist(fr), unlist(to))), n = sum(Wl), directed = FALSE)
  list(g = g, vid = vid, W = Wl, w = unlist(wt))
}

## Dijkstra from one OR SEVERAL start cells: one row per start cell, one column per passable cell (map units; Inf = unreachable)
graph_dist_matrix <- function(gg, rows, cols) {
  s <- gg$vid[cbind(rows, cols)]; stopifnot(all(s > 0))
  igraph::distances(gg$g, v = s, weights = gg$w)
}
## distance to the NEAREST start cell, as a matrix shaped like the grid; unreachable -> NA
graph_dist_field <- function(gg, row, col) {
  m <- graph_dist_matrix(gg, row, col)
  d <- if (nrow(m) == 1) m[1, ] else Reduce(pmin, lapply(seq_len(nrow(m)), function(i) m[i, ]))
  d[is.infinite(d)] <- NA
  out <- matrix(NA_real_, nrow(gg$W), ncol(gg$W)); out[gg$W] <- d
  out
}

## A start REGION: a set of cells (two-column matrix, row and column in the grid) that all count as distance 0. One virtual vertex is joined
## to each of them by a zero-length edge, so ONE Dijkstra search from it gives the distance to the nearest cell of the region. `gg` itself is
## left as it is - the extended copy belongs to this region only - so a region can never act as a shortcut for any other start.
graph_region_source <- function(gg, cells) {
  .wd_need("igraph")
  stopifnot(is.matrix(cells), ncol(cells) == 2, nrow(cells) >= 1)
  v <- gg$vid[cells]; stopifnot(all(v > 0))
  n <- igraph::vcount(gg$g)
  list(g = igraph::add_edges(igraph::add_vertices(gg$g, 1), as.vector(rbind(n + 1L, v))),
       w = c(gg$w, rep(0, length(v))), src = n + 1L, n = n)
}
## distance (map units) from the region to every real vertex of the grid graph; Inf = unreachable
graph_region_dist <- function(rg) igraph::distances(rg$g, v = rg$src, weights = rg$w)[1, seq_len(rg$n)]

## ---- 2. snap points to the nearest passable cell, within a tolerance ---------------------------------------
snap_to_passable <- function(xy, r, W, tol_m) {
  nr <- nrow(W); nc <- ncol(W); res <- terra::res(r)[1]; k <- ceiling(tol_m / res)
  out <- data.frame(row = NA_integer_, col = NA_integer_, snap_m = NA_real_)[rep(1, nrow(xy)), ]
  rc <- terra::rowColFromCell(r, terra::cellFromXY(r, xy))
  for (i in seq_len(nrow(xy))) {
    if (anyNA(rc[i, ])) next
    rows <- max(1, rc[i, 1] - k):min(nr, rc[i, 1] + k); cols <- max(1, rc[i, 2] - k):min(nc, rc[i, 2] + k)
    sub <- which(W[rows, cols, drop = FALSE], arr.ind = TRUE)
    if (!nrow(sub)) next
    d <- sqrt((terra::xFromCol(r, cols[sub[, 2]]) - xy[i, 1])^2 + (terra::yFromRow(r, rows[sub[, 1]]) - xy[i, 2])^2)
    j <- which.min(d)
    if (d[j] <= tol_m) out[i, ] <- list(rows[sub[j, 1]], cols[sub[j, 2]], d[j])
  }
  out
}

## ---- 2b. start boxes ---------------------------------------------------------------------------------------
## One box c(xmin, xmax, ymin, ymax) = lon_min, lon_max, lat_min, lat_max (WGS84 degrees), or a four-column matrix with one row per box
## -> numeric matrix with columns xmin, xmax, ymin, ymax. NULL -> a matrix with no rows. `what` names the setting in the error message.
.wd_box_matrix <- function(b, what = "`ref_box`") {
  cols <- c("xmin", "xmax", "ymin", "ymax")
  if (is.null(b)) return(matrix(numeric(0), 0, 4, dimnames = list(NULL, cols)))
  bad <- function() stop(what, " must be a box c(xmin, xmax, ymin, ymax) = lon_min, lon_max, lat_min, lat_max in WGS84 degrees ",
                         "(lon within -180..180, lat within -90..90, min < max), or a matrix with one such row per box; got ",
                         paste(format(unlist(b)), collapse = ", "), call. = FALSE)
  if ((is.matrix(b) || is.data.frame(b)) && ncol(b) != 4) bad()
  if (!(is.matrix(b) || is.data.frame(b)) && length(b) %% 4 != 0) bad()
  m <- if (is.matrix(b) || is.data.frame(b)) matrix(as.numeric(as.matrix(b)), ncol = 4) else matrix(as.numeric(b), ncol = 4, byrow = TRUE)
  if (!nrow(m) || anyNA(m) || any(!is.finite(m)) || any(abs(m[, 1:2]) > 180) || any(abs(m[, 3:4]) > 90) ||
      any(m[, 1] >= m[, 2]) || any(m[, 3] >= m[, 4])) bad()
  dimnames(m) <- list(NULL, cols)
  m
}

## The outline of a box as a closed ring of lon/lat vertices, a vertex about every 0.05 degree (3-6 km): the straight pieces between the
## vertices then follow the parallels and meridians through any map projection (the deviation is metres), and geodesics between them
## stay on the parallel.
.wd_box_ring <- function(b, step = 0.05) {
  nx <- min(3000, max(2, ceiling((b[2] - b[1]) / step) + 1)); ny <- min(3000, max(2, ceiling((b[4] - b[3]) / step) + 1))
  seg <- function(x1, y1, x2, y2, n) cbind(seq(x1, x2, length.out = n), seq(y1, y2, length.out = n))
  rbind(seg(b[1], b[3], b[2], b[3], nx), seg(b[2], b[3], b[2], b[4], ny), seg(b[2], b[4], b[1], b[4], nx), seg(b[1], b[4], b[1], b[3], ny))
}

## The boxes as a SpatVector of polygons (WGS84), one per row of the matrix
.wd_box_vect <- function(boxm) {
  polys <- lapply(seq_len(nrow(boxm)), function(i) terra::vect(.wd_box_ring(boxm[i, ]), type = "polygons", crs = "EPSG:4326"))
  if (length(polys) == 1) polys[[1]] else do.call(rbind, polys)
}

## The passable cells that ONE box touches, on the grid `r` (the box as a projected polygon, `W` = passable matrix). A cell counts when the
## box passes through any part of it, not only its centre: a box thinner than a cell (the HELCOM line is 220 m) would otherwise often
## miss every cell centre. Returns list(cells = two-column matrix of row and column of the passable cells, touched = cells touched in all).
.wd_box_cells <- function(poly_m, r, W) {
  none <- list(cells = matrix(integer(0), 0, 2), touched = 0L)
  res <- terra::res(r)[1]; e <- terra::ext(poly_m); er <- terra::ext(r)
  x1 <- max(er[1], e[1] - 2 * res); x2 <- min(er[2], e[2] + 2 * res); y1 <- max(er[3], e[3] - 2 * res); y2 <- min(er[4], e[4] + 2 * res)
  if (x1 >= x2 || y1 >= y2) return(none)                           # the box lies outside the grid
  win <- terra::crop(terra::rast(r), terra::ext(x1, x2, y1, y2))   # same grid, only the window around the box: cheap on a large grid
  hit <- which(!is.na(terra::values(terra::rasterize(poly_m, win, field = 1, touches = TRUE))[, 1]))
  if (!length(hit)) return(none)
  rc <- terra::rowColFromCell(r, terra::cellFromXY(r, terra::xyFromCell(win, hit)))
  list(cells = rc[W[rc], , drop = FALSE], touched = nrow(rc))
}

## ---- 3. the main function ----------------------------------------------------------------------------------
## points   : data.frame / matrix, columns lon, lat (WGS84). One row per UNIQUE sampling location is enough.
## ref      : the start POINT(s): c(lon, lat) for one, or a two-column matrix (one row per point) for several, e.g. all the
##            entrances to a sea. NULL for none (then `ref_box` is needed). Each start point must lie within snap_tol_km of a passable cell.
## ref_box  : the start BOX(es): c(xmin, xmax, ymin, ymax) = lon_min, lon_max, lat_min, lat_max in WGS84 degrees for one box, or a
##            four-column matrix with one row per box. NULL for none (then `ref` is needed). A box is a region: every passable cell it
##            touches is a start cell at distance 0, so a site's route distance is the shortest route to ANY of them, and a site inside the
##            box has distance 0. Boxes may be thinner than a cell (the HELCOM line is 220 m high): they still select the cells they pass
##            through. The route ends at the CENTRE of such a cell, so a box is resolved to about half a cell (0.5 km at 1 km cells).
##            Cells of the box that are not passable (land, for domain ocean) are not used; a box that touches no passable cell stops
##            with an error.
##            With several starts of either kind the distance is the shortest route to the NEAREST one; table$ref_id says which: the
##            start points in the order given, then the boxes numbered on from there.
## polygons : SpatVector of OCEAN polygons (e.g. the Natural Earth `ocean` you already load). Any CRS is accepted
##            (converted to lon/lat and cropped internally), but pass the lon/lat original if you still have it.
## domain   : "ocean" = route only through water; "land" = route only over land.
## polygons_domain: what `polygons` delineate: "ocean" (default) or "land". With the default, domain "land" routes over the part OUTSIDE the
##            polygons (inverted mask). Pass "land" together with domain "land" when `polygons` are land polygons (the Natural Earth
##            layer the wingen Rmd loads for domain_type "land"): the route then runs inside them.
## crs      : "auto" = local Lambert azimuthal equal-area centred on the data (distortion < ~1 % within ~1500 km),
##            or any metric CRS such as "EPSG:3035". Never use lon/lat or Mollweide here.
## min_cover: a cell is passable only if >= this fraction of it is inside the domain. 0.99 = "any land touching
##            the cell makes it land": thin spits / peninsulas stay barriers, but straits narrower than one cell
##            close (they then show up as status "unreachable", not as silent errors). 0.5 = majority rule:
##            keeps narrow straits open but can leak across thin land.
## pad_km   : margin of the grid around the points and the starts. A route that would have to leave the grid cannot be
##            found, so if a site is "unreachable" while its water body continues past the grid edge (an enclosed lagoon
##            does not), the margin is doubled and the search repeated, up to max_pad_km. settings$pad_km /
##            settings$pad_attempts report what was used.
## snap_tol_km: a point further than this from a passable cell gets NA + status "no_passable_cell_near".
## paths    : TRUE additionally returns the route of every point as a SpatVector of lines (for QA plots).
## negative_side: optional, off by default. Sites on the "other side" of the start (e.g. Kattegat/Skagerrak when the start is the Oresund)
##            get NEGATIVE signed distances. Give a lon/lat box c(xmin, xmax, ymin, ymax) or a polygon SpatVector. The result then
##            has a column signed_route_km (= route_km, negated inside the box/polygon). Without it there is no such column and
##            every distance is >= 0: the right choice when the study area lies on one side of the start, as inside the Baltic Sea.
## Returns list(table, field, paths, settings, box_cells). One row per input point, in input order:
##   route_km    the shortest route through the passable cells (water or land, see `domain`) from the nearest start: THE METRIC
##   geodesic_km the AS-THE-CROW-FLIES distance (geodesic on the ellipsoid) from the site to that same start, ignoring land: to the start
##               point, or to the nearest edge of the start box (0 inside it). For sanity checks
##   circuity    route_km / geodesic_km: should not be below 1 (a route cannot be shorter than the crow-fly line; ~1 % grid noise aside, and
##               a cell or two for sites within a few cells of the start). Large is not an error - it says "check this route"; compare
##               with neighbouring sites. NA when geodesic_km is 0
##   snap_km     how far the coordinate was moved to the centre of the passable cell used as its node (0 - ~0.7 cell if its own cell is passable)
##   status      "ok"; "no_passable_cell_near" = no passable cell within snap_tol_km (coordinate on land in this coastline, or its water
##               body is missing from it); "unreachable" = snapped fine, but cut off from every start (landlocked, or a passage narrower than a cell)
##   signed_route_km  only with `negative_side`: route_km, made NEGATIVE for points inside it; the value to use as a predictor when the data
##               lie on both sides of the start
##   ref_id      which start the route_km was measured from: 1..n for the start points (rows of `ref`), then n+1.. for the boxes
## field = the distance raster (km) on the local grid; paths = route lines (paths = TRUE); settings = what was actually used, including
## reference_box_info (one row per box: its edges, the cells it touches and how many of them are passable);
## box_cells = for each box, the lon/lat of the centres of the passable cells it stands for (a list of two-column matrices).
water_route_km <- function(points, ref = NULL, polygons, domain = c("ocean", "land"), res_km = 1,
                           crs = "auto", min_cover = 0.99, connectivity = 16, pad_km = 200, max_pad_km = 800,
                           snap_tol_km = max(3 * res_km, 3), paths = FALSE, negative_side = NULL,
                           ref_box = NULL, polygons_domain = "ocean") {
  domain <- match.arg(domain)
  stopifnot(length(polygons_domain) == 1, polygons_domain %in% c("ocean", "land"))
  run <- function(pad) .water_route_once(points, ref, polygons, domain = domain, res_km = res_km, crs = crs,
                                         min_cover = min_cover, connectivity = connectivity, pad_km = pad,
                                         snap_tol_km = snap_tol_km, paths = paths, negative_side = negative_side,
                                         ref_box = ref_box, polygons_domain = polygons_domain)
  pad <- pad_km; attempts <- 1; res <- run(pad)
  while (length(res$suspect) && pad * 2 <= max_pad_km) {
    bigger <- try(run(pad * 2), silent = TRUE)          # too many cells -> keep what we have
    if (inherits(bigger, "try-error")) break
    pad <- pad * 2; attempts <- attempts + 1; res <- bigger
  }
  if (length(res$suspect))
    warning(length(res$suspect), " site(s) are cut off from the start although their water body continues past the grid ",
            "edge (margin ", pad, " km): the route may need a larger margin than max_pad_km.", call. = FALSE)
  res$settings$pad_km <- pad; res$settings$pad_attempts <- attempts
  res$suspect <- NULL
  res
}

## One search on a grid with the given margin. Also returns $suspect: unreachable sites whose connected water body touches
## the grid edge, i.e. sites where the margin may be the reason they are unreachable.
.wd_ref_matrix <- function(ref) {
  if (is.null(ref)) return(matrix(numeric(0), 0, 2))
  if (!is.matrix(ref) && !is.data.frame(ref) && length(ref) %% 2 != 0) stop("`ref` must be c(lon, lat) or a two-column matrix of lon/lat", call. = FALSE)
  m <- if (is.matrix(ref) || is.data.frame(ref)) matrix(as.numeric(as.matrix(ref)), ncol = 2) else matrix(as.numeric(ref), ncol = 2, byrow = TRUE)
  if (!nrow(m) || anyNA(m) || any(abs(m[, 1]) > 180) || any(abs(m[, 2]) > 90)) stop("`ref` must be WGS84 lon/lat degrees (lon in [-180, 180], lat in [-90, 90])", call. = FALSE)
  m
}

.water_route_once <- function(points, ref = NULL, polygons, domain = c("ocean", "land"), res_km = 1,
                              crs = "auto", min_cover = 0.99, connectivity = 16, pad_km = 200,
                              snap_tol_km = max(3 * res_km, 3), paths = FALSE, negative_side = NULL,
                              ref_box = NULL, polygons_domain = "ocean") {
  .wd_need()
  domain <- match.arg(domain)
  pts <- as.matrix(points[, 1:2]); colnames(pts) <- c("lon", "lat")
  refm <- .wd_ref_matrix(ref); npt <- nrow(refm)
  boxm <- .wd_box_matrix(ref_box); nbx <- nrow(boxm); nref <- npt + nbx
  if (nref == 0) stop("No start: give start point(s) in `ref` and/or start box(es) in `ref_box`.", call. = FALSE)
  stopifnot(is.numeric(pts), all(is.finite(pts)), all(abs(pts[, 1]) <= 180), all(abs(pts[, 2]) <= 90), inherits(polygons, "SpatVector"))
  corners <- if (nbx) do.call(rbind, lapply(seq_len(nbx), function(i) cbind(boxm[i, c(1, 2, 2, 1)], boxm[i, c(3, 3, 4, 4)]))) else matrix(numeric(0), 0, 2)
  allxy <- rbind(refm, unname(corners), pts)
  if (diff(range(allxy[, 1])) > 180)
    stop("Longitude range > 180 degrees (global data or data crossing the antimeridian): not supported by this single-grid method.")
  if (identical(crs, "auto")) {
    crs <- sprintf("+proj=laea +lat_0=%.4f +lon_0=%.4f +datum=WGS84 +units=m +no_defs",
                   mean(range(allxy[, 2])), mean(range(allxy[, 1])))
  }
  res_m <- res_km * 1000
  xy  <- terra::crds(terra::project(terra::vect(pts, crs = "EPSG:4326"), crs))
  rxy <- if (npt) terra::crds(terra::project(terra::vect(refm, crs = "EPSG:4326"), crs)) else matrix(numeric(0), 0, 2)
  box_ll <- if (nbx) .wd_box_vect(boxm) else NULL
  box_m  <- if (nbx) terra::project(box_ll, crs) else NULL
  bxy    <- if (nbx) terra::crds(box_m) else matrix(numeric(0), 0, 2)
  pad <- pad_km * 1000; allm <- rbind(unname(rxy), unname(bxy), unname(xy))
  span_km <- sqrt(diff(range(allm[, 1]))^2 + diff(range(allm[, 2]))^2) / 1000
  if (span_km > 3500) warning("Study area spans ~", round(span_km), " km: one metric grid distorts distances by several %. ",
                              "Split the data into regions or pass a regional `crs`.")
  e <- terra::ext(floor((min(allm[, 1]) - pad) / res_m) * res_m, ceiling((max(allm[, 1]) + pad) / res_m) * res_m,
                  floor((min(allm[, 2]) - pad) / res_m) * res_m, ceiling((max(allm[, 2]) + pad) / res_m) * res_m)
  r <- terra::rast(e, resolution = res_m, crs = crs)
  if (terra::ncell(r) > 6e7) stop("Grid would have ", terra::ncell(r), " cells: increase res_km or reduce pad_km / study area.")

  # fraction of each cell inside the POLYGONS -> passable matrix (inverted when the polygons delineate the other medium)
  pl <- if (terra::is.lonlat(polygons)) polygons else terra::project(polygons, "EPSG:4326")
  dlat <- 1.5 * pad_km / 111 + 1                                   # crop in lon/lat FIRST: projecting a global polygon is slow/unsafe
  dlon <- 1.5 * pad_km / (111 * cos(min(85, max(abs(allxy[, 2])) + dlat) * pi / 180)) + 1
  o <- terra::crop(pl, terra::ext(max(-180, min(allxy[, 1]) - dlon), min(180, max(allxy[, 1]) + dlon),
                                  max(-90, min(allxy[, 2]) - dlat), min(90, max(allxy[, 2]) + dlat)))
  o <- terra::crop(terra::project(o, crs), e)
  frac <- terra::values(terra::rasterize(o, r, cover = TRUE))[, 1]; frac[is.na(frac)] <- 0
  if (domain != polygons_domain) frac <- 1 - frac
  W <- matrix(frac >= min_cover, nrow(r), ncol(r), byrow = TRUE)
  if (!any(W)) stop("No passable cells: check `polygons` (ocean polygons expected), `domain` and `crs`.")

  sn <- snap_to_passable(rbind(rxy, xy), r, W, snap_tol_km * 1000)
  ref_sn <- sn[seq_len(npt), , drop = FALSE]; s <- sn[npt + seq_len(nrow(xy)), , drop = FALSE]
  bad_ref <- which(is.na(ref_sn$row))
  if (length(bad_ref))
    stop("Start point(s) ", paste(sprintf("#%d (lon %.4f, lat %.4f)", bad_ref, refm[bad_ref, 1], refm[bad_ref, 2]), collapse = ", "),
         " lie more than ", snap_tol_km, " km from a passable cell. A narrow passage can be missing from the grid at this cell size ",
         "(e.g. the Little Belt at 1 km): use a finer res_km, move the point into open water, or raise snap_tol_km.", call. = FALSE)

  # the passable cells each start box touches
  seeds <- vector("list", nbx); touched <- integer(nbx)
  for (j in seq_len(nbx)) {
    sj <- .wd_box_cells(box_m[j, ], r, W); touched[j] <- sj$touched; seeds[[j]] <- sj$cells
    if (!nrow(sj$cells))
      stop("Start box #", npt + j, " (lon ", format(boxm[j, 1]), " to ", format(boxm[j, 2]), ", lat ", format(boxm[j, 3]), " to ", format(boxm[j, 4]),
           ") touches no passable cell: ", if (touched[j]) paste0("the ", touched[j], " grid cell(s) it touches are not passable (", if (domain == "ocean") "land" else "water", ", or too little of the cell lies in the ", domain, ")")
                                           else "it lies outside the grid or the coastline that was loaded (the coastline must cover the box)",
           ". Check the coordinates (lon first, then lat) or use a finer res_km.", call. = FALSE)
  }

  gg <- grid_graph(W, res_m, connectivity, corner_safe = TRUE)
  dl <- vector("list", nref)                                        # one distance vector (km) per start, over the passable cells
  if (npt) { mp <- graph_dist_matrix(gg, ref_sn$row, ref_sn$col) / 1000; for (i in seq_len(npt)) dl[[i]] <- mp[i, ]; rm(mp) }
  for (j in seq_len(nbx)) dl[[npt + j]] <- graph_region_dist(graph_region_source(gg, seeds[[j]])) / 1000
  dmin <- Reduce(pmin, dl)                                          # distance to the NEAREST start
  dmin[is.infinite(dmin)] <- NA
  D <- matrix(NA_real_, nrow(W), ncol(W)); D[gg$W] <- dmin
  hit <- !is.na(s$row)
  route <- rep(NA_real_, nrow(s)); route[hit] <- D[cbind(s$row[hit], s$col[hit])]
  ref_id <- rep(NA_integer_, nrow(s))
  if (any(hit)) {
    vh <- gg$vid[cbind(s$row[hit], s$col[hit])]
    ms <- matrix(vapply(dl, function(d) d[vh], numeric(length(vh))), ncol = nref)     # one row per hit site, one column per start
    nearest <- max.col(-ms, ties.method = "first")
    nearest[!is.finite(apply(ms, 1, min))] <- NA_integer_            # not reachable from any start
    ref_id[hit] <- nearest
  }
  pv <- terra::vect(pts, crs = "EPSG:4326")
  dg <- matrix(NA_real_, nrow(pts), nref)                            # crow-flies distance (km) of every site to every start
  if (npt) dg[, seq_len(npt)] <- as.matrix(terra::distance(pv, terra::vect(refm, crs = "EPSG:4326"))) / 1000
  if (nbx) dg[, npt + seq_len(nbx)] <- as.matrix(terra::distance(pv, box_ll)) / 1000
  geodesic <- dg[cbind(seq_len(nrow(dg)), ifelse(is.na(ref_id), max.col(-dg, ties.method = "first"), ref_id))]
  tab <- data.frame(lon = pts[, 1], lat = pts[, 2], route_km = route, geodesic_km = geodesic,
                    circuity = ifelse(geodesic > 0, route / geodesic, NA_real_), snap_km = s$snap_m / 1000,
                    status = ifelse(is.na(s$row), "no_passable_cell_near", ifelse(is.na(route), "unreachable", "ok")))
  if (!is.null(negative_side)) {
    neg <- if (inherits(negative_side, "SpatVector")) {
      rowSums(terra::relate(terra::vect(pts, crs = "EPSG:4326"), terra::project(negative_side, "EPSG:4326"), "within")) > 0
    } else {
      stopifnot(is.numeric(negative_side), length(negative_side) == 4)
      pts[, 1] >= negative_side[1] & pts[, 1] <= negative_side[2] & pts[, 2] >= negative_side[3] & pts[, 2] <= negative_side[4]
    }
    tab$signed_route_km <- ifelse(neg, -tab$route_km, tab$route_km)
  }
  tab$ref_id <- ref_id
  field <- terra::setValues(terra::rast(r), as.vector(t(D)))

  suspect <- integer(0)                       # unreachable, but their water body continues past the grid edge?
  unreach <- which(hit & is.na(route))
  if (length(unreach)) {
    comp <- igraph::components(gg$g)$membership
    ring <- matrix(FALSE, nrow(W), ncol(W)); ring[1, ] <- TRUE; ring[nrow(W), ] <- TRUE; ring[, 1] <- TRUE; ring[, ncol(W)] <- TRUE
    touching <- unique(comp[gg$vid[ring & gg$W]])
    suspect <- unreach[comp[gg$vid[cbind(s$row[unreach], s$col[unreach])]] %in% touching]
  }

  route_lines <- NULL
  if (paths && any(tab$status == "ok")) {          # nothing routable -> no lines (an empty set cannot be built)
    ok <- which(tab$status == "ok")
    tv <- gg$vid[cbind(s$row[ok], s$col[ok])]; src_v <- if (npt) gg$vid[cbind(ref_sn$row, ref_sn$col)] else integer(0)
    vp <- vector("list", length(ok))
    for (k in unique(ref_id[ok])) { sel <- which(ref_id[ok] == k)             # each site's route goes to ITS nearest start
      if (k <= npt) {
        vp[sel] <- igraph::shortest_paths(gg$g, from = src_v[k], to = tv[sel], weights = gg$w, output = "vpath")$vpath
      } else {                                                                  # from the box: the virtual first vertex is not a cell
        rg <- graph_region_source(gg, seeds[[k - npt]])
        vp[sel] <- lapply(igraph::shortest_paths(rg$g, from = rg$src, to = tv[sel], weights = rg$w, output = "vpath")$vpath,
                          function(p) as.integer(p)[-1])
      }
    }
    lin <- which(gg$W); nr <- nrow(W)
    coords <- lapply(vp, function(p) { l <- lin[as.integer(p)]
      if (length(l) == 1) l <- c(l, l)                                          # a site in its start cell: a line of length 0
      cbind(terra::xFromCol(r, (l - 1) %/% nr + 1), terra::yFromRow(r, (l - 1) %% nr + 1)) })
    route_lines <- terra::vect(coords, type = "lines", crs = crs); route_lines$point_row <- ok
  }

  to_lonlat <- function(rows, cols) unname(terra::crds(terra::project(terra::vect(cbind(terra::xFromCol(r, cols), terra::yFromRow(r, rows)), crs = crs), "EPSG:4326")))
  box_cells <- lapply(seeds, function(cells) to_lonlat(cells[, 1], cells[, 2]))
  box_info <- NULL
  if (nbx) box_info <- data.frame(box = npt + seq_len(nbx), xmin = boxm[, 1], xmax = boxm[, 2], ymin = boxm[, 3], ymax = boxm[, 4],
                                  cells_touched = touched, cells_passable = vapply(seeds, nrow, integer(1)),
                                  passable_lon_min = vapply(box_cells, function(m) min(m[, 1]), numeric(1)),
                                  passable_lon_max = vapply(box_cells, function(m) max(m[, 1]), numeric(1)),
                                  passable_lat_min = vapply(box_cells, function(m) min(m[, 2]), numeric(1)),
                                  passable_lat_max = vapply(box_cells, function(m) max(m[, 2]), numeric(1)))
  list(table = tab, field = field, paths = route_lines, suspect = suspect, box_cells = box_cells,
       settings = list(domain = domain, polygons_domain = polygons_domain, res_km = res_km, crs = crs, min_cover = min_cover,
                       connectivity = connectivity, pad_km = pad_km, snap_tol_km = snap_tol_km, passable_cells = sum(W),
                       edges = igraph::ecount(gg$g), n_references = nref, n_reference_points = npt, n_reference_boxes = nbx,
                       reference_snapped_km = ref_sn$snap_m / 1000,
                       reference_snapped_lonlat = if (npt) to_lonlat(ref_sn$row, ref_sn$col) else matrix(numeric(0), 0, 2),
                       reference_box_info = box_info,
                       negative_side = if (is.null(negative_side)) "none" else if (inherits(negative_side, "SpatVector")) "polygon" else as.numeric(negative_side)))
}

## ---- 4. helpers for the wingen pipeline (used by code/2_Herring_Wingen.Rmd, section 2.7) --------------------

## One row per unique sampling site. The key is paste(Latitude, Longitude): the same key the Rmd's
## assign_distance_clusters() uses to map individuals to sites.
water_dist_sites <- function(coords_df, lon = "Longitude", lat = "Latitude") {
  key <- paste(coords_df[[lat]], coords_df[[lon]])
  keep <- !duplicated(key)
  data.frame(site_key = key[keep], lon = as.numeric(coords_df[[lon]][keep]), lat = as.numeric(coords_df[[lat]][keep]),
             stringsAsFactors = FALSE)
}

## One or several boxes -> numeric matrix, one row per box, columns xmin, xmax, ymin, ymax (WGS84 degrees); empty / "none" -> NULL.
## Boxes are separated by semicolons. A box is written with two labelled intervals, in either order, in the notation of the study design
##     "lon [10.60, 12.00], lat [57.7395, 57.7415]"      (also "lon 10.60 12.00, lat 57.7395 57.7415"; = longitude 10.60-12.00 E, latitude 57.7395-57.7415 N)
## so that latitude and longitude cannot be mixed up, or as the four bare numbers "xmin,xmax,ymin,ymax" = lon_min, lon_max, lat_min, lat_max.
## `what` names the setting in the error message.
water_dist_parse_boxes <- function(x, what = "WATER_DIST_REFERENCE_BOXES") {
  x <- trimws(x)
  if (!length(x) || is.na(x) || !nzchar(x) || tolower(x) %in% c("none", "null", "na")) return(NULL)
  fail <- function(part) stop(what, " must list boxes separated by semicolons, each written like \"lon [10.60, 12.00], lat [57.7395, 57.7415]\" ",
                              "or as xmin,xmax,ymin,ymax (= lon_min, lon_max, lat_min, lat_max) in WGS84 degrees with min < max; got '", part, "'", call. = FALSE)
  num <- "[-+]?(?:[0-9]+(?:\\.[0-9]+)?|\\.[0-9]+)"
  lon_word <- "\\b(?:longitude|long|lng|lon)\\b"; lat_word <- "\\b(?:latitude|lat)\\b"
  one <- function(p) {
    v <- as.numeric(regmatches(p, gregexpr(num, p, perl = TRUE))[[1]])
    rest <- gsub(num, "", p, perl = TRUE)
    rest <- gsub(paste0("(?i)", lon_word, "|", lat_word), "", rest, perl = TRUE)
    if (length(v) != 4 || anyNA(v) || nzchar(gsub("[][[:space:],=:()]", "", rest))) fail(p)
    lon_at <- regexpr(lon_word, p, perl = TRUE, ignore.case = TRUE); lat_at <- regexpr(lat_word, p, perl = TRUE, ignore.case = TRUE)
    if (lon_at < 0 && lat_at < 0) return(v)                         # four bare numbers: xmin, xmax, ymin, ymax
    if (lon_at < 0 || lat_at < 0) fail(p)                            # one label without the other
    if (lon_at < lat_at) v else v[c(3, 4, 1, 2)]                     # the interval written first belongs to the label written first
  }
  parts <- trimws(strsplit(x, ";")[[1]]); parts <- parts[nzchar(parts)]
  if (!length(parts)) fail(x)
  .wd_box_matrix(do.call(rbind, lapply(parts, one)), what)
}

## Exactly one box -> numeric(4) c(xmin, xmax, ymin, ymax); empty / "none" -> NULL (used for the optional negative side)
water_dist_parse_box <- function(x) {
  m <- water_dist_parse_boxes(x, what = "WATER_DIST_NEGATIVE_BOX")
  if (is.null(m)) return(NULL)
  if (nrow(m) != 1) stop("WATER_DIST_NEGATIVE_BOX takes exactly one box; got ", nrow(m), call. = FALSE)
  unname(m[1, ])
}

## One-off QA figure for checking that routes stay in the medium: the distance field with 200 km contours, the route to each
## site (BLACK where the signed distance is positive, RED where it is negative, i.e. on the far side of the entrance - only when a negative
## side is set), the start point(s) as gold triangles and the start box(es) as gold lines, the negative-side box (dashed red) when one is
## set, and a purple cross on every site that could not
## be routed. At most max_routes routes are drawn: the farthest, the most winding (circuity) and the furthest-snapped sites always,
## the rest evenly thinned - deterministically, with NO random numbers, so the caller's RNG state (e.g. the neutral simulations)
## is untouched. Drawn on the local grid the routes were computed on. Needs a res from water_route_km(..., paths = TRUE).
water_dist_qa_plot <- function(res, ref, file, max_routes = 300, ref_box = NULL) {
  .wd_need("terra")
  fld <- res$field; tab <- res$table; crs_ <- terra::crs(fld); refm <- .wd_ref_matrix(ref); boxm <- .wd_box_matrix(ref_box)
  npt <- nrow(refm); nbx <- nrow(boxm)
  ex <- terra::ext(fld); px_w <- 1500; px_h <- max(700, min(2400, round(px_w * (ex[4] - ex[3]) / (ex[2] - ex[1]))))
  tmp <- paste0(file, ".tmp", Sys.getpid(), ".png")
  grDevices::png(tmp, width = px_w, height = px_h, res = 170)
  dev <- grDevices::dev.cur(); closed <- FALSE
  on.exit({ if (!closed && dev %in% grDevices::dev.list()) grDevices::dev.off(dev); unlink(tmp) }, add = TRUE)
  ok_rows <- which(tab$status == "ok")
  signed <- "signed_route_km" %in% names(tab)
  neg <- if (signed) !is.na(tab$signed_route_km) & tab$signed_route_km < 0 else rep(FALSE, nrow(tab))
  red <- "#e31a1c"
  starts <- c(if (npt) sprintf("%d start point%s", npt, if (npt > 1) "s" else ""), if (nbx) sprintf("%d start box%s", nbx, if (nbx > 1) "es" else ""))
  # The figure is opened with terra::plot of the distance field (blocked cells beige), not with an empty graphics::plot frame that terra
  # then draws into: terra keeps one clip rectangle per graphics device, and only a terra::plot that does not add sets it for the
  # device in use. Drawing into a bare frame would clip this figure to whatever an earlier terra plot on the same device slot
  # (an earlier figure, a knitr chunk) left behind.
  terra::plot(fld, col = grDevices::hcl.colors(60, "YlGnBu", rev = TRUE), colNA = "#e9e4d8", legend = FALSE, axes = FALSE, box = FALSE,
              mar = c(1, 1, 3, 1), maxcell = 4e6,
              main = sprintf("Shortest %s routes to the nearest of %s: %d sites routed%s, %d not  [%g km cells]",
                             res$settings$domain, paste(starts, collapse = " and "), length(ok_rows),
                             if (signed) sprintf(" (%d on the negative side)", sum(neg & tab$status == "ok")) else "",
                             nrow(tab) - length(ok_rows), res$settings$res_km), cex.main = 0.7)
  rng <- range(terra::values(fld), na.rm = TRUE)
  levs <- seq(200, floor(max(abs(rng)) / 200) * 200, by = 200)
  if (length(levs)) try(terra::contour(fld, levels = levs, add = TRUE, col = "#00000055", lwd = 0.6, labcex = 0.7), silent = TRUE)
  ns <- res$settings$negative_side
  if (is.numeric(ns) && length(ns) == 4) {                                               # the negative-side box, densified so it bends with the projection
    seg <- function(x1, y1, x2, y2) cbind(seq(x1, x2, length.out = 60), seq(y1, y2, length.out = 60))
    ring <- rbind(seg(ns[1], ns[3], ns[2], ns[3]), seg(ns[2], ns[3], ns[2], ns[4]), seg(ns[2], ns[4], ns[1], ns[4]), seg(ns[1], ns[4], ns[1], ns[3]))
    try(terra::lines(terra::project(terra::vect(list(ring), type = "lines", crs = "EPSG:4326"), crs_), col = red, lty = 2, lwd = 1.3), silent = TRUE)
  }
  if (!is.null(res$paths) && nrow(res$paths) > 0) {
    idx <- seq_len(nrow(res$paths)); rows <- res$paths$point_row
    if (length(idx) > max_routes) {
      top <- function(x, n) order(x, decreasing = TRUE)[seq_len(min(n, length(x)))]
      always <- unique(c(top(tab$route_km[rows], 25), top(tab$circuity[rows], 25), top(tab$snap_km[rows], 25)))
      rest <- setdiff(idx, always)
      idx <- sort(unique(c(always, rest[round(seq(1, length(rest), length.out = max_routes - length(always)))])))
    }
    is_neg <- neg[rows[idx]]
    terra::lines(res$paths[idx, ], col = "#ffffffaa", lwd = 2.4)                         # white halo: readable on light and dark water
    if (any(!is_neg)) terra::lines(res$paths[idx[!is_neg], ], col = "#1a1a1a", lwd = 0.9)
    if (any(is_neg))  terra::lines(res$paths[idx[is_neg], ],  col = red, lwd = 1.2)
  }
  proj <- function(df) terra::project(terra::vect(as.matrix(df[, c("lon", "lat")]), crs = "EPSG:4326"), crs_)
  pos_ok <- ok_rows[!neg[ok_rows]]; neg_ok <- ok_rows[neg[ok_rows]]
  if (length(pos_ok)) terra::points(proj(tab[pos_ok, ]), pch = 21, bg = "white", col = "#1a1a1a", cex = 0.55, lwd = 0.7)
  if (length(neg_ok)) terra::points(proj(tab[neg_ok, ]), pch = 21, bg = red, col = "#7f0000", cex = 0.6, lwd = 0.7)
  bad <- which(tab$status != "ok")
  if (length(bad)) terra::points(proj(tab[bad, ]), pch = 4, col = "#7b3294", cex = 1.2, lwd = 1.8)
  gold_box <- "#e6b800"
  for (j in seq_len(nbx)) {                                                              # a start box is a line at this scale: black casing, gold core
    bl <- terra::project(terra::vect(list(.wd_box_ring(boxm[j, ])), type = "lines", crs = "EPSG:4326"), crs_)
    terra::lines(bl, col = "black", lwd = 4.2); terra::lines(bl, col = gold_box, lwd = 2.2)
    if (npt + nbx > 1) { cc <- terra::crds(terra::project(terra::vect(cbind(mean(boxm[j, 1:2]), boxm[j, 4]), crs = "EPSG:4326"), crs_))
                         graphics::text(cc[1], cc[2], labels = npt + j, pos = 3, cex = 0.8, font = 2, offset = 0.7) }
  }
  if (npt) {
    rp <- terra::project(terra::vect(refm, crs = "EPSG:4326"), crs_)
    terra::points(rp, pch = 24, bg = "gold", col = "black", cex = 1.6, lwd = 1.2)
    if (npt + nbx > 1) graphics::text(terra::crds(rp)[, 1], terra::crds(rp)[, 2], labels = seq_len(npt), pos = 3, cex = 0.8, font = 2, offset = 0.7)
  }
  L <- list(); add <- function(label, pch = NA_real_, lty = NA_real_, lwd = 1, col = "black", bg = NA_character_) L[[length(L) + 1]] <<- list(label, pch, lty, lwd, col, bg)
  if (npt) add(if (npt + nbx > 1) "start points (numbered = ref_id)" else "start point", pch = 24, bg = "gold")
  if (nbx) add(if (npt + nbx > 1) "start boxes (numbered = ref_id; distance 0)" else "start box (distance 0)", lty = 1, lwd = 3, col = gold_box)
  add(if (signed) "site, positive side" else "site", pch = 21, col = "#1a1a1a", bg = "white")
  add(if (signed) "route, positive side" else "route", lty = 1, col = "#1a1a1a")
  if (length(neg_ok) || any(neg)) { add("site / route, negative side", pch = 21, lty = 1, col = red, bg = red); add("negative-side box", lty = 2, col = red) }
  if (length(bad)) add("site NOT routed", pch = 4, col = "#7b3294")
  add(sprintf("contours every 200 km (%g-%g km)", round(rng[1]), round(rng[2])), lty = 1, col = "#00000055")
  pick <- function(i) vapply(L, function(a) a[[i]], if (i == 1) "" else if (i %in% c(5, 6)) "" else numeric(1))
  graphics::legend(x = ex[1] + 0.015 * (ex[2] - ex[1]), y = ex[4] - 0.012 * (ex[4] - ex[3]), xjust = 0, yjust = 1, xpd = NA,
                   legend = pick(1), pch = pick(2), lty = pick(3), lwd = pick(4), col = pick(5), pt.bg = pick(6), bg = "#ffffffcc", box.col = NA, cex = 0.72)
  grDevices::dev.off(dev); closed <- TRUE
  file.rename(tmp, file)
  invisible(file)
}

## "lon,lat; lon,lat; ..." (WGS84 degrees) -> numeric matrix with one row per start point; empty -> NULL
water_dist_parse_points <- function(x) {
  x <- trimws(x)
  if (!nzchar(x) || tolower(x) %in% c("none", "null", "na")) return(NULL)
  parts <- trimws(strsplit(x, ";")[[1]]); parts <- parts[nzchar(parts)]
  m <- t(vapply(parts, function(p) { v <- suppressWarnings(as.numeric(strsplit(p, "[, ]+")[[1]]))
                                     if (length(v) == 2 && !anyNA(v)) v else c(NA_real_, NA_real_) }, numeric(2)))
  if (!length(parts) || anyNA(m) || any(abs(m[, 1]) > 180) || any(abs(m[, 2]) > 90))
    stop("WATER_DIST_REFERENCE_POINTS must be 'lon,lat' pairs in WGS84 degrees separated by semicolons, ",
         "e.g. \"12.80,55.75; 10.95,55.45\"; got '", x, "'", call. = FALSE)
  m <- unname(m); dimnames(m) <- NULL
  m
}

## Compute-or-load. The cache file name carries a hash of EVERY input that can change the numbers (start points and boxes, grid
## settings, coastline fingerprint, the set of sites, algorithm version), so a changed setting can never reuse a stale file.
## Written atomically (temp file + rename): genes that start at the same moment on a cold cache each compute it once
## (seconds) and none ever reads a half-written file.
## Returns list(table = one row per site incl. site_key, meta, file, from_cache).
water_dist_cached <- function(sites, ref = NULL, polygons, cache_dir, domain = "ocean", res_km = 1, min_cover = 0.99,
                              snap_tol_km = max(3 * res_km, 3), negative_side = NULL, pad_km = 200,
                              connectivity = 16, coastline_id = "", qa_plot = FALSE, ref_box = NULL, polygons_domain = "ocean") {
  .wd_need(c("terra", "igraph", "jsonlite"))
  settings <- list(version = WATER_DIST_VERSION, ref = as.numeric(.wd_ref_matrix(ref)),
                   ref_box = if (is.null(ref_box)) NA else as.numeric(.wd_box_matrix(ref_box)),
                   domain = domain, polygons_domain = polygons_domain, res_km = res_km,
                   min_cover = min_cover, snap_tol_km = snap_tol_km,
                   negative_side = if (is.null(negative_side)) NA else as.numeric(negative_side),
                   pad_km = pad_km, connectivity = connectivity, coastline = coastline_id, sites = sites$site_key)
  key_file <- tempfile(fileext = ".txt"); on.exit(unlink(key_file), add = TRUE)
  writeLines(paste(deparse(settings, control = c("keepNA", "digits17")), collapse = "\n"), key_file)
  hash <- substr(unname(tools::md5sum(key_file)), 1, 12)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  tsv <- file.path(cache_dir, sprintf("water_distance_%s.tsv", hash)); meta_file <- sub("\\.tsv$", ".json", tsv)
  qa_file <- sub("\\.tsv$", "_routes.png", tsv)
  have_cache <- file.exists(tsv) && file.exists(meta_file)
  need_plot <- isTRUE(qa_plot) && !file.exists(qa_file)        # the figure is made once, next to the cache file
  read_cached <- function() {
    meta <- jsonlite::read_json(meta_file, simplifyVector = TRUE)
    meta$reference_snapped_lonlat <- matrix(as.numeric(unlist(meta$reference_snapped_lonlat)), ncol = 2)   # 0 rows when there is no start point
    meta$reference_snapped_km <- as.numeric(unlist(meta$reference_snapped_km))
    if (!is.data.frame(meta$reference_box_info)) meta$reference_box_info <- NULL                             # no start box
    list(table = utils::read.delim(tsv, stringsAsFactors = FALSE, quote = "", na.strings = "NA",
                                   colClasses = c(site_key = "character", status = "character")),
         meta = meta)
  }
  if (have_cache && !need_plot)
    return(c(read_cached(), list(file = tsv, from_cache = TRUE, qa_file = if (file.exists(qa_file)) qa_file else NA_character_)))
  res <- water_route_km(sites[, c("lon", "lat")], ref, polygons, domain = domain, res_km = res_km, min_cover = min_cover,
                        connectivity = connectivity, pad_km = pad_km, snap_tol_km = snap_tol_km, negative_side = negative_side,
                        paths = need_plot, ref_box = ref_box, polygons_domain = polygons_domain)
  if (have_cache) {
    out <- read_cached()                                         # keep what is on disk; only the figure is new
  } else {
    tab <- data.frame(site_key = sites$site_key, res$table, stringsAsFactors = FALSE)
    columns <- list(                                              # the sidecar explains every column of the TSV
      site_key = "paste(Latitude, Longitude), the key the Rmd uses to join sites to individuals",
      lon = "longitude of the site (WGS84)", lat = "latitude of the site (WGS84)",
      route_km = "THE METRIC: shortest route from the nearest start (a start point, or any passable cell of a start box), through passable cells only (water for domain ocean, land for domain land), in km; always >= 0, and 0 for a site inside a start box",
      geodesic_km = "AS-THE-CROW-FLIES distance (geodesic on the ellipsoid, ignores land) from the site to the start route_km was measured from: the start point, or the nearest edge of the start box (0 inside it)",
      circuity = "route_km / geodesic_km; should not be below 1 (about 1 % grid noise aside, and a cell or two close to the start); NA when geodesic_km is 0; a large value is not an error, just worth a look",
      snap_km = "how far the coordinate was moved to the centre of the passable grid cell used as its node (includes the offset inside its own cell)",
      status = "ok | no_passable_cell_near (no passable cell within the snap tolerance) | unreachable (placed, but cut off from every start)")
    if (!is.null(negative_side))
      columns$signed_route_km <- "route_km, made NEGATIVE for sites inside WATER_DIST_NEGATIVE_BOX (the far side of the start); this is then the predictor"
    columns$ref_id <- "which start was nearest by route: the row number in WATER_DIST_REFERENCE_POINTS, then the boxes of WATER_DIST_REFERENCE_BOXES numbered on from the last point"
    meta <- c(res$settings, list(hash = hash, version = WATER_DIST_VERSION, coastline = coastline_id,
                                 n_sites = nrow(sites), created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), columns = columns))
    pid <- Sys.getpid(); tsv_tmp <- paste0(tsv, ".tmp", pid); meta_tmp <- paste0(meta_file, ".tmp", pid)
    utils::write.table(tab, tsv_tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    jsonlite::write_json(meta, meta_tmp, auto_unbox = TRUE, digits = NA, na = "null")
    file.rename(meta_tmp, meta_file); file.rename(tsv_tmp, tsv)   # TSV last: its presence means the metadata exists too
    out <- list(table = tab, meta = meta)
  }
  if (need_plot) {                                                # a plotting problem must never break a run
    made <- try(water_dist_qa_plot(res, ref, qa_file, ref_box = ref_box), silent = TRUE)
    if (inherits(made, "try-error")) warning("The water-route QA figure could not be written: ", conditionMessage(attr(made, "condition")), call. = FALSE)
  }
  c(out, list(file = tsv, from_cache = have_cache, qa_file = if (file.exists(qa_file)) qa_file else NA_character_))
}

## NULL if every site was routed; otherwise a readable list of the sites that were not (for stop() / message()).
water_dist_failure_text <- function(tab, snap_tol_km = NA, n_show = 15) {
  bad <- tab[tab$status != "ok", , drop = FALSE]
  if (!nrow(bad)) return(NULL)
  shown <- utils::head(bad, n_show)
  what <- ifelse(is.na(shown$snap_km),
                 if (is.na(snap_tol_km)) "no passable cell within the snap tolerance" else sprintf("no passable cell within %.1f km", snap_tol_km),
                 sprintf("snapped by %.1f km, but cut off from every start", shown$snap_km))
  lines <- sprintf("  %-21s (lon %.5f, lat %.5f)  %s", shown$status, shown$lon, shown$lat, what)
  paste0(nrow(bad), " of ", nrow(tab), " sampling site(s) got no water-route distance:\n", paste(lines, collapse = "\n"),
         if (nrow(bad) > n_show) paste0("\n  ... and ", nrow(bad) - n_show, " more") else "",
         "\nno_passable_cell_near = no passable cell within the snap tolerance (coordinate on land in this coastline, or a water body missing from it);",
         "\nunreachable = snapped, but cut off from every start point and box (e.g. a strait narrower than one grid cell). ",
         "Raise WATER_DIST_SNAP_TOL_KM, use a finer WATER_DIST_RES_KM or a finer coastline; or set WATER_DIST_ON_FAILED_SITES: drop.")
}

## Per-cluster summary of a per-individual predictor (named by Sample ID). The mean is the SAME rule as the Rmd's
## env_for_cluster(): mean over the member individuals. range = max - min over the members: because two sites' route
## distances cannot differ by more than the water distance between them, a range much larger than the cluster's own
## straight-line diameter means the (straight-line) cluster straddles a land barrier. flag_km = threshold for that.
water_dist_cluster_summary <- function(cluster_groups, ind_values, flag_km) {
  rows <- lapply(cluster_groups, function(ids) {
    v <- ind_values[ids]; ok <- is.finite(v)
    if (!any(ok)) return(c(mean = NA_real_, min = NA_real_, max = NA_real_, range = NA_real_, n_missing = length(ids)))
    c(mean = mean(v[ok]), min = min(v[ok]), max = max(v[ok]), range = max(v[ok]) - min(v[ok]), n_missing = sum(!ok))
  })
  out <- as.data.frame(do.call(rbind, rows), row.names = NULL)
  out$flag <- !is.na(out$range) & out$range > flag_km
  out
}
