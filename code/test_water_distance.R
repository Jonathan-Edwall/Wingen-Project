# test_water_distance.R  -  run:  Rscript test_water_distance.R   (set NE_OCEAN_GEOJSON=/path/ne_10m_ocean.geojson
# to also run the real-coastline checks; Parts A, H, J and K need no data at all)
suppressPackageStartupMessages(library(terra))   # the tests use terra directly; water_distance.R itself attaches nothing
source("water_distance.R")
n_checks <- 0L
check <- function(label, ok) { n_checks <<- n_checks + 1L; cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label)); if (!isTRUE(ok)) quit(status = 1) }
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
err_text <- function(expr) { r <- try(expr, silent = TRUE); if (inherits(r, "try-error")) conditionMessage(attr(r, "condition")) else "" }

cat("Part A - toy grids (no data needed)\n")
# A1. a 1-cell-thick DIAGONAL land wall must separate the two sides
n <- 9; W <- matrix(TRUE, n, n); diag(W) <- FALSE
d <- graph_dist_field(grid_graph(W, 1000, 16, corner_safe = TRUE), 1, n)
check("diagonal land wall is NOT crossed (other side unreachable)", is.na(d[n, 1]))
d_leaky <- graph_dist_field(grid_graph(W, 1000, 8, corner_safe = FALSE), 1, n)   # = what terra::gridDist does
check("without the corner rule the same wall leaks (reproduces gridDist)", is.finite(d_leaky[n, 1]))
# A2. open-water accuracy of the 16-neighbour graph (square cells)
N <- 801; c0 <- 401; Wopen <- matrix(TRUE, N, N); D <- graph_dist_field(grid_graph(Wopen, 1, 16), c0, c0)
ang <- seq(0, 45, by = 2.5) * pi / 180
err <- sapply(ang, function(a) { r_ <- c0 - round(300 * sin(a)); c_ <- c0 + round(300 * cos(a)); D[r_, c_] / sqrt((r_ - c0)^2 + (c_ - c0)^2) - 1 })
check(sprintf("16-neighbour open-water error: max +%.1f %% (< 3 %%), mean +%.1f %%", 100 * max(err), 100 * mean(err)), max(err) < 0.03 && min(err) > -1e-9)
# A3. a start REGION: a set of cells at distance 0, searched once from a virtual vertex
n <- 41; Wr <- matrix(TRUE, n, n); gg <- grid_graph(Wr, 1000, 16); nv <- igraph::vcount(gg$g); ne <- igraph::ecount(gg$g)
rg <- graph_region_source(gg, cbind(1:n, 1L))                                   # the whole western column is the start
Dr <- matrix(NA_real_, n, n); Dr[gg$W] <- graph_region_dist(rg)
check("start region: distance 0 on every cell of the region, 20 km twenty cells away, and the grid graph itself is left untouched",
      all(Dr[, 1] == 0) && all(abs(Dr[, 21] - 20000) < 1e-9) && igraph::vcount(gg$g) == nv && igraph::ecount(gg$g) == ne)
# A4. two separate water bodies (a land wall across row 21): the region joins them, a point start does not
Ww <- matrix(TRUE, n, n); Ww[21, ] <- FALSE; gw <- grid_graph(Ww, 1000, 16)
rw <- graph_region_source(gw, rbind(c(1L, 1L), c(41L, 1L)))                      # one start cell on each side of the wall
Dw <- matrix(NA_real_, n, n); Dw[gw$W] <- graph_region_dist(rw)
Dp <- graph_dist_field(gw, 1, 1)                                                 # a point start in the north
check("start region: both sides of the wall reached from their own start cell; a point start in the north does not leak into the south",
      is.finite(Dw[40, 2]) && is.finite(Dw[2, 2]) && Dw[40, 2] < 5000 && is.na(Dp[40, 2]) && is.finite(Dp[2, 2]))

cat("Part H - start boxes on a synthetic coastline (no data needed)\n")
bx <- function(x) unname(water_dist_parse_boxes(x)[1, ])
helcom <- c(10.60, 12.00, 57.7395, 57.7415)
check("water_dist_parse_boxes: labelled in either order (the notation of the study design), bare numbers, several boxes, empty",
      identical(bx("lon [10.60, 12.00], lat [57.7395, 57.7415]"), helcom) && identical(bx("lat [57.7395, 57.7415], lon [10.60, 12.00]"), helcom) &&
      identical(bx("latitude = [57.7395, 57.7415] longitude = [10.60, 12.00]"), helcom) && identical(bx("lon 10.60 12.00, lat 57.7395 57.7415"), helcom) &&
      identical(bx("10.60,12.00,57.7395,57.7415"), helcom) && identical(bx("  LON [10.60, 12.00],  LAT [57.7395, 57.7415] ; "), helcom))
two <- water_dist_parse_boxes("lon [10.60, 12.00], lat [57.7395, 57.7415]; lon [-5.5, -3], lat [48, 49.5]")
check("water_dist_parse_boxes: two boxes give two rows (negative longitudes allowed); empty and 'none' give NULL",
      identical(dim(two), c(2L, 4L)) && isTRUE(all.equal(unname(two[2, ]), c(-5.5, -3, 48, 49.5))) &&
      is.null(water_dist_parse_boxes("")) && is.null(water_dist_parse_boxes("none")))
bad_boxes <- c("lon [12, 10], lat [57, 58]", "lon [10, 12]", "lon [10, 12], lon [3, 4]", "lat [100, 101], lon [1, 2]", "10,12,57",
               "lon [10..12], lat [1..2]", "lon [10, 12], lat [57, 57]", "lon [190, 200], lat [1, 2]", "somewhere")
check("water_dist_parse_boxes: reversed or empty interval, a missing or doubled label, latitude beyond 90, three numbers, '..', unknown words are refused",
      all(vapply(bad_boxes, function(b) fails(water_dist_parse_boxes(b)), logical(1))))
check("water_dist_parse_box (the optional negative side): one box only", identical(water_dist_parse_box("lon [10, 12], lat [57, 58]"), c(10, 12, 57, 58)) &&
      fails(water_dist_parse_box("lon [1, 2], lat [3, 4]; lon [5, 6], lat [7, 8]")))
check("the message names the setting and the accepted notation", grepl("WATER_DIST_REFERENCE_BOXES", err_text(water_dist_parse_boxes("lon [10, 12]"))) &&
      grepl("lon \\[10.60, 12.00\\], lat \\[57.7395, 57.7415\\]", err_text(water_dist_parse_boxes("lon [10, 12]"))))

# a 10 x 6 degree sea with a square island (4-6 E, 52-54 N); the start box is a 220 m high strip along 55 N from 2 E to 8 E, north of the island
sea <- vect("POLYGON ((0 50, 10 50, 10 56, 0 56, 0 50), (4 52, 6 52, 6 54, 4 54, 4 52))", crs = "EPSG:4326")
box <- c(2, 8, 55.0, 55.002)
sites <- data.frame(lon = c(3, 5, 5, 8, 5, 5.0), lat = c(54.5, 51, 55.001, 52, 55.5, 53.0))   # west of the island, south of it, IN the box, south-east, north of the box, ON the island
o <- water_route_km(sites, ref = NULL, polygons = sea, ref_box = box, res_km = 1, paths = TRUE)
tb <- o$table
way <- function(p) sum(vapply(seq_len(nrow(p) - 1), function(i) as.numeric(distance(vect(p[i, , drop = FALSE], crs = "EPSG:4326"), vect(p[i + 1, , drop = FALSE], crs = "EPSG:4326"))), numeric(1))) / 1000
detour <- way(rbind(c(5, 51), c(4, 52), c(4, 54), c(4, 55)))          # hug the island's west side, then straight north to the strip
check(sprintf("a start box needs no start point (ref = NULL); the route of a site is the shortest water route to the box: %.1f km (crow flies %.1f km)", tb$route_km[1], tb$geodesic_km[1]),
      all(tb$status[1:5] == "ok") && abs(tb$route_km[1] - tb$geodesic_km[1]) < 1.5 && abs(tb$route_km[5] - tb$geodesic_km[5]) < 1.5)
check(sprintf("the island is a barrier: south of it the route is %.0f km, within 2 %% of the path around its west corner (%.0f km), well above the crow-fly line (%.0f km)",
              tb$route_km[2], detour, tb$geodesic_km[2]),
      abs(tb$route_km[2] / detour - 1) < 0.02 && tb$route_km[2] > tb$geodesic_km[2] * 1.03)
check("a site inside the box has distance 0 (and no circuity: 0 / 0); a site on the island is flagged, not routed",
      tb$route_km[3] == 0 && tb$geodesic_km[3] == 0 && is.na(tb$circuity[3]) && tb$status[6] == "no_passable_cell_near")
check("a box is reached from the north as well as from the south: no sign, both sides positive", all(tb$route_km[tb$status == "ok"] >= 0) && !("signed_route_km" %in% names(tb)))
ring <- .wd_box_ring(box, step = 0.01)
d_edge <- apply(as.matrix(distance(vect(as.matrix(sites[1:5, ]), crs = "EPSG:4326"), vect(ring, crs = "EPSG:4326"))), 1, min) / 1000
check("geodesic_km of a box = distance to its nearest edge (0 inside), agrees with a brute-force search along the outline (< 0.05 km)",
      all(abs(tb$geodesic_km[c(1, 2, 4, 5)] - d_edge[c(1, 2, 4, 5)]) < 0.05))
bi <- o$settings$reference_box_info
check(sprintf("settings$reference_box_info: edges, %d cells touched, all passable; box_cells holds their lon/lat; the distance field is 0 on all of them and positive 3 cells away",
              bi$cells_touched),
      nrow(bi) == 1 && bi$cells_touched == bi$cells_passable && nrow(o$box_cells[[1]]) == bi$cells_passable && bi$box == 1 &&
      all(extract(o$field, project(vect(o$box_cells[[1]], crs = "EPSG:4326"), crs(o$field)))[, 2] == 0) &&
      extract(o$field, project(vect(cbind(5, 55.05), crs = "EPSG:4326"), crs(o$field)))[1, 2] > 3)
check("route polylines (also for the site in the box: a line of length 0) have the reported length (+-0.1 %)",
      all(abs(perim(o$paths) / 1000 - tb$route_km[o$paths$point_row]) < 1e-3 * pmax(tb$route_km[o$paths$point_row], 1)))
# ---- start points and a box together: ref_id counts the points first, then the boxes
m2 <- water_route_km(rbind(sites[c(1, 4, 5), ], data.frame(lon = 9.5, lat = 50.5)), ref = c(9.5, 51.0), polygons = sea, ref_box = box, res_km = 1)
check("a point and a box together: the nearest of the two wins; ref_id 1 = the point, 2 = the first box; geodesic_km refers to the one that won",
      identical(m2$table$ref_id, c(2L, 1L, 2L, 1L)) && m2$table$geodesic_km[4] < 80 && m2$table$geodesic_km[2] > 140 && m2$table$geodesic_km[2] < 160 && m2$settings$n_references == 2 &&
      m2$settings$n_reference_points == 1 && m2$settings$n_reference_boxes == 1 && identical(m2$settings$reference_box_info$box, 2L))
# ---- a box is the same as a start point at every one of its cells (an independent route to the same numbers); a coarse grid keeps the number of starts small
smallsea <- vect("POLYGON ((0 50, 6 50, 6 54, 0 54, 0 50))", crs = "EPSG:4326")
bsmall <- c(1, 5, 52.0, 52.003); ss <- data.frame(lon = c(0.6, 3, 5.4, 2.2, 3.3), lat = c(50.4, 50.3, 53.5, 52.0005, 53.9))
laea <- "+proj=laea +lat_0=52 +lon_0=3 +datum=WGS84 +units=m +no_defs"
ob <- water_route_km(ss, polygons = smallsea, ref_box = bsmall, res_km = 4, crs = laea, pad_km = 150)
op <- water_route_km(ss, ref = ob$box_cells[[1]], polygons = smallsea, res_km = 4, crs = laea, pad_km = 150)
check(sprintf("box = one start point on the centre of each of its %d cells: identical distances to the nearest cell (largest difference %.1e km)", nrow(ob$box_cells[[1]]),
              max(abs(ob$table$route_km - op$table$route_km))),
      all(ob$table$status == "ok") && max(abs(ob$table$route_km - op$table$route_km)) < 1e-6)
# ---- a thin box resolves to about half a cell wherever it lies relative to the cell centres (a cell-centre rule would often miss a 100 m strip)
lat0 <- seq(52.0, 52.0 + 0.0100, length.out = 12); south <- data.frame(lon = 3, lat = 50.5); dev <- numeric(length(lat0))
for (i in seq_along(lat0)) {
  oi <- water_route_km(south, polygons = smallsea, ref_box = c(1, 5, lat0[i], lat0[i] + 0.001), res_km = 2, crs = laea, pad_km = 100)
  dev[i] <- oi$table$route_km - oi$table$geodesic_km
}
check(sprintf("a 110 m high box moved in 12 steps across one cell row: the distance always stays within one 2 km cell of the crow-fly line (largest deviation %.2f km)", max(abs(dev))),
      all(abs(dev) < 2.5))
check("a box on land (the island, as a start) and a box outside the coastline stop with a message saying why and how to fix it",
      grepl("touches no passable cell", err_text(water_route_km(sites, polygons = sea, ref_box = c(4.3, 5.7, 52.5, 53.5)))) &&
      grepl("outside the grid|touches no passable", err_text(water_route_km(sites, polygons = sea, ref_box = c(20, 21, 60, 61)))))
check("neither a start point nor a start box stops with a message", grepl("No start", err_text(water_route_km(sites, polygons = sea))))
check("a start box with the wrong shape or order stops (xmin >= xmax, three numbers)",
      fails(water_route_km(sites, polygons = sea, ref_box = c(8, 2, 55, 55.002))) && fails(water_route_km(sites, polygons = sea, ref_box = c(2, 8, 55))))
cat("Part H2 - caching with a start box\n")
cs <- water_dist_sites(data.frame(`Sample ID` = paste0("s", 1:5), Latitude = ss$lat, Longitude = ss$lon, check.names = FALSE))
cd <- tempfile("wdbox"); on.exit(unlink(cd, recursive = TRUE), add = TRUE)
c1 <- water_dist_cached(cs, polygons = smallsea, cache_dir = cd, res_km = 4, ref_box = water_dist_parse_boxes("lat [52, 52.003], lon [1, 5]"), pad_km = 150)
c2 <- water_dist_cached(cs, polygons = smallsea, cache_dir = cd, res_km = 4, ref_box = bsmall, pad_km = 150)
c3 <- water_dist_cached(cs, polygons = smallsea, cache_dir = cd, res_km = 4, ref_box = c(1, 5, 52.0, 52.004), pad_km = 150)
check("cache with a box only: computed once, then loaded from disk with the same table, whichever way the box was written; another box gets its own file",
      !c1$from_cache && c2$from_cache && c1$file == c2$file && isTRUE(all.equal(c1$table$route_km, c2$table$route_km)) && !c3$from_cache && c3$file != c1$file)
check("cache: no signed_route_km column or description unless a negative side is set; ref_id and the description of every column are there",
      !("signed_route_km" %in% names(c1$table)) && !("signed_route_km" %in% names(c2$meta$columns)) && all(c("route_km", "geodesic_km", "circuity", "ref_id") %in% names(c2$meta$columns)) &&
      identical(names(c2$table), c("site_key", "lon", "lat", "route_km", "geodesic_km", "circuity", "snap_km", "status", "ref_id")))
check("cache: the sidecar of a box-only run reads back with 0 start points, 1 box and its table of cells; a run with start points reads back the snapped positions",
      c2$meta$n_reference_points == 0 && c2$meta$n_reference_boxes == 1 && is.data.frame(c2$meta$reference_box_info) && c2$meta$reference_box_info$cells_passable > 50 &&
      identical(dim(c2$meta$reference_snapped_lonlat), c(0L, 2L)) &&
      identical(dim(water_dist_cached(cs, ref = c(3, 50.5), polygons = smallsea, cache_dir = cd, res_km = 4, pad_km = 150)$meta$reference_snapped_lonlat), c(1L, 2L)))
qd <- tempfile("wdboxqa"); on.exit(unlink(qd, recursive = TRUE), add = TRUE)
q <- water_dist_cached(cs, polygons = smallsea, cache_dir = qd, res_km = 4, ref_box = bsmall, pad_km = 150, qa_plot = TRUE)
check("QA figure with a start box and no start point: written as a real PNG", !is.na(q$qa_file) && file.exists(q$qa_file) && file.size(q$qa_file) > 15000)
invisible(file.copy(q$qa_file, file.path(Sys.getenv("WD_TEST_FIGURES_DIR", tempdir()), "box_only_qa_figure.png"), overwrite = TRUE))

cat("Part H3 - a route figure is not clipped to an earlier plot\n")
# terra keeps one clip rectangle per graphics device. The figure above has just been drawn on the first device slot, so a second figure with
# a much larger extent must still show its whole distance field (in a bare graphics frame it was cut to the extent of the first one).
wide_sea <- vect(rbind(c(-5, 45), c(40, 45), c(40, 70), c(-5, 70), c(-5, 45)), type = "polygons", crs = "EPSG:4326")
wide_sites <- water_dist_sites(data.frame(`Sample ID` = c("a", "b", "c"), Latitude = c(58, 60, 61), Longitude = c(20, 25, 18), check.names = FALSE))
field_share <- function(png_file) {   # share of the figure frame (beige blocked cells + field colours) that shows the distance field
  v <- suppressWarnings(as.matrix(terra::values(terra::rast(png_file))))[, 1:3, drop = FALSE]
  beige <- abs(v[, 1] - 233) < 3 & abs(v[, 2] - 228) < 3 & abs(v[, 3] - 216) < 3
  white <- v[, 1] > 250 & v[, 2] > 250 & v[, 3] > 250
  sum(!beige & !white) / sum(!white)
}
qw <- water_dist_cached(wide_sites, polygons = wide_sea, cache_dir = tempfile("wdwide"), res_km = 8, ref_box = c(14, 16, 55.4, 55.402), qa_plot = TRUE)
check(sprintf("a second QA figure on the same device slot shows its whole field (%.0f %% of its frame; %.0f %% for the small sea of the first figure)",
              100 * field_share(qw$qa_file), 100 * field_share(q$qa_file)),
      !is.na(qw$qa_file) && field_share(qw$qa_file) > 0.8 && field_share(q$qa_file) < 0.5)

cat("Part J - the signed distance is optional\n")
sg <- water_route_km(sites, polygons = sea, ref_box = box, res_km = 1, negative_side = c(0, 10, 55.0015, 56))   # everything north of the box's centre line
okr <- sg$table$status == "ok"
check("with negative_side the table gets signed_route_km = route_km, negated inside the side (here the site north of the box); without it there is no such column",
      "signed_route_km" %in% names(sg$table) && all(abs(abs(sg$table$signed_route_km[okr]) - sg$table$route_km[okr]) < 1e-9) &&
      sg$table$signed_route_km[5] < 0 && all(sg$table$signed_route_km[okr & seq_along(okr) != 5] >= 0) &&
      !("signed_route_km" %in% names(water_route_km(sites, polygons = sea, ref_box = box, res_km = 1)$table)) &&
      isTRUE(all.equal(sg$table$route_km, tb$route_km)))

cat("Part K - land domain with land polygons\n")
island <- vect("POLYGON ((4 52, 6 52, 6 54, 4 54, 4 52))", crs = "EPSG:4326")                      # the island itself, as LAND polygons
sl <- data.frame(lon = c(4.5, 5.5, 5.0, 3.0), lat = c(52.5, 53.5, 53.0, 53.0))                      # three sites on the island, one in the sea
a <- water_route_km(sl, c(4.2, 52.2), island, domain = "land", polygons_domain = "land", res_km = 1, pad_km = 80)$table
b <- water_route_km(sl, c(4.2, 52.2), sea, domain = "land", res_km = 1, pad_km = 80)$table              # the same island, as the hole of the OCEAN polygons
check("domain land: land polygons (polygons_domain = 'land') and ocean polygons (the default) give the same land routes; the site in the sea is flagged",
      identical(a$status, c("ok", "ok", "ok", "no_passable_cell_near")) && isTRUE(all.equal(a$route_km, b$route_km)) && identical(b$status, a$status) && a$route_km[2] > a$geodesic_km[2] * 0.99)
a2 <- water_route_km(data.frame(lon = c(5, 8), lat = c(51, 55)), c(9, 51), island, domain = "ocean", polygons_domain = "land", res_km = 1, pad_km = 80)$table
b2 <- water_route_km(data.frame(lon = c(5, 8), lat = c(51, 55)), c(9, 51), sea, domain = "ocean", res_km = 1, pad_km = 80)$table
check("domain ocean with land polygons (route outside them) equals domain ocean with ocean polygons", a2$status[1] == "ok" && isTRUE(all.equal(a2$route_km[1], b2$route_km[1])))
check("an unknown polygons_domain is refused", fails(water_route_km(sl, c(4.2, 52.2), island, domain = "land", polygons_domain = "sea")))

cat("Part B - real coastline (Natural Earth 10 m ocean)\n")
f <- Sys.getenv("NE_OCEAN_GEOJSON", "")
if (!nzchar(f) || !file.exists(f)) { cat("  skipped: set NE_OCEAN_GEOJSON\n"); cat(sprintf("all %d checks passed (synthetic parts only)\n", n_checks)); quit(status = 0) }
ocean <- vect(f); if (is.na(crs(ocean)) || crs(ocean) == "") crs(ocean) <- "EPSG:4326"     # deliberately NOT pre-cropped
st <- read.csv(text = "name,lon,lat
Kattegat,11.5,57.0
GreatBelt,10.95,55.35
GdanskBay,19.0,54.7
PuckBay_outer,18.60,54.80
PuckBay_inner,18.55,54.68
VistulaLagoon,19.9,54.55
CuronianLagoon,21.1,55.2
BothnianBay,23.0,65.0")
ref <- c(12.80, 55.75)                                                                   # Oresund
out <- water_route_km(st[, c("lon", "lat")], ref, ocean, res_km = 1, paths = TRUE,       # crs = "auto"
                      negative_side = c(8, 13.0, 55.95, 60))                              # Kattegat/Skagerrak side of Oresund
tb <- out$table; rownames(tb) <- st$name; print(round(tb[, c("route_km", "geodesic_km", "circuity", "snap_km")], 1)); print(tb$status)
check("Hel spit stays a barrier: inner Puck Bay is > 30 km farther than the outer coast",
      tb["PuckBay_inner", "route_km"] - tb["PuckBay_outer", "route_km"] > 30)
check("Vistula Lagoon (not in the ocean polygon, 7 km from open water) is flagged, not silently snapped",
      tb["VistulaLagoon", "status"] == "no_passable_cell_near")
check("Curonian Lagoon (strait narrower than one 1-km cell) is reported 'unreachable'", tb["CuronianLagoon", "status"] == "unreachable")
ok <- tb$status == "ok"
check("every route is at least as long as the straight line", all(tb$route_km[ok] >= tb$geodesic_km[ok] * 0.995))
pl_km <- perim(out$paths) / 1000
check("extracted route polylines have the same length as the reported distances (+-0.1 %)",
      all(abs(pl_km / tb$route_km[tb$status == "ok"] - 1) < 1e-3))
check("negative_side: Kattegat gets a negative signed distance, Baltic sites stay positive",
      tb["Kattegat", "signed_route_km"] < 0 && all(tb[setdiff(rownames(tb), "Kattegat"), "signed_route_km"][ok[-1]] > 0))
out3035 <- water_route_km(st[, c("lon", "lat")], ref, ocean, res_km = 1, crs = "EPSG:3035")$table
check("local LAEA ('auto') and EPSG:3035 agree within 1 %", all(abs(out3035$route_km / tb$route_km - 1) < 0.01, na.rm = TRUE))

cat("Part C - same machinery, terrestrial domain (land-only route)\n")
land <- water_route_km(data.frame(lon = 18.07, lat = 59.33), ref = c(13.0, 55.6), polygons = ocean, domain = "land", res_km = 1)
cat(sprintf("  Malmo -> Stockholm over land: %.0f km (straight %.0f km)\n", land$table$route_km, land$table$geodesic_km))
check("land route exists and is 1.0-1.3 x the straight line", land$table$status == "ok" && land$table$circuity > 1 && land$table$circuity < 1.3)
cat("Part D - outside the Baltic ('auto' local projection) and guard rails\n")
ad <- water_route_km(data.frame(lon = 13.7, lat = 45.6), c(18.5, 40.0), ocean)$table          # Otranto -> Trieste
check(sprintf("Adriatic route %.0f km vs %.0f km straight", ad$route_km, ad$geodesic_km), ad$status == "ok" && ad$circuity > 1 && ad$circuity < 1.3)
hw <- water_route_km(data.frame(lon = -154.9, lat = 19.9), c(-157.9, 21.3), ocean)$table       # Honolulu -> off Hilo
check(sprintf("Hawaii route %.0f km vs %.0f km straight", hw$route_km, hw$geodesic_km), hw$status == "ok" && hw$circuity > 1 && hw$circuity < 1.3)
check("data crossing the antimeridian is refused, not silently wrong",
      inherits(try(water_route_km(data.frame(lon = -170, lat = 10), c(170, 10), ocean), silent = TRUE), "try-error"))
cat("Part E - pipeline helpers (cache, cluster summary, failure report)\n")
coords <- data.frame(`Sample ID` = paste0("ind", 1:8), check.names = FALSE,
                     Latitude  = c(54.80, 54.80, 54.68, 54.68, 55.35, 55.35, 54.55, 65.0),
                     Longitude = c(18.60, 18.60, 18.55, 18.55, 10.95, 10.95, 19.9, 23.0))   # outer Puck x2, inner Puck x2, Great Belt x2, Vistula Lagoon, Bothnian Bay
sites <- water_dist_sites(coords)
check("water_dist_sites: one row per unique coordinate", nrow(sites) == 5 && !anyDuplicated(sites$site_key))
check("water_dist_parse_box", identical(water_dist_parse_box("8, 13, 55.95, 60"), c(8, 13, 55.95, 60)) && is.null(water_dist_parse_box("")) &&
      inherits(try(water_dist_parse_box("13,8,1,2"), silent = TRUE), "try-error"))
cdir <- tempfile("wdcache"); on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
a <- water_dist_cached(sites, ref, ocean, cdir, coastline_id = "ne-test")
b <- water_dist_cached(sites, ref, ocean, cdir, coastline_id = "ne-test")
check("cache: first call computes, second call loads the same table from disk", !a$from_cache && b$from_cache &&
      isTRUE(all.equal(a$table$route_km, b$table$route_km)) && identical(a$table$site_key, b$table$site_key))
check("cache: atomic write leaves exactly one .tsv + one .json and no temp files", setequal(sub("^water_distance_[0-9a-f]{12}", "", list.files(cdir)), c(".tsv", ".json")))
d <- water_dist_cached(sites, ref, ocean, cdir, res_km = 2, coastline_id = "ne-test")
check("cache: a changed setting gets its own file (never a stale hit)", !d$from_cache && d$file != a$file)
ids_by_cluster <- list(`1` = c("ind1", "ind2", "ind3", "ind4"),   # outer + inner Puck Bay: ~14 km apart in a straight line, ~50 km by water
                       `2` = c("ind5", "ind6"),                    # same site twice
                       `3` = "ind7", `4` = "ind8")
ind_val <- setNames(a$table$route_km[match(paste(coords$Latitude, coords$Longitude), a$table$site_key)], coords$`Sample ID`)
cs <- water_dist_cluster_summary(ids_by_cluster, ind_val, flag_km = 15)
check(sprintf("cluster summary flags the Hel-spit cluster (range %.0f km) but not a one-site cluster", cs$range[1]),
      isTRUE(cs$flag[1]) && !cs$flag[2] && cs$range[2] == 0)
check("cluster summary: all-NA cluster gives NA mean, no flag", is.na(cs$mean[3]) && !cs$flag[3] && cs$n_missing[3] == 1)
txt <- water_dist_failure_text(a$table, snap_tol_km = 3)
check("failure report names the lagoon site and the tolerance", !is.null(txt) && grepl("no_passable_cell_near", txt) && grepl("1 of 5", txt) && grepl("within 3.0 km", txt))
check("failure report is NULL when everything routed", is.null(water_dist_failure_text(a$table[a$table$status == "ok", ])))
check("water_distance.R attaches nothing (no masking of Rmd functions)", !any(c("package:igraph") %in% search()))
small <- sites[1:3, ]                                              # the QA-figure checks need only a few sites
qdir <- tempfile("wdqa"); on.exit(unlink(qdir, recursive = TRUE), add = TRUE)
q1 <- water_dist_cached(small, ref, ocean, qdir, qa_plot = TRUE, coastline_id = "ne-test")
check("QA figure: written once, next to the cache file, as a real PNG", !is.na(q1$qa_file) && file.exists(q1$qa_file) && file.size(q1$qa_file) > 20000)
m1 <- file.mtime(q1$qa_file); Sys.sleep(1.1)
q2 <- water_dist_cached(small, ref, ocean, qdir, qa_plot = TRUE, coastline_id = "ne-test")
check("QA figure: a cache hit neither recomputes nor redraws, and leaves no temp files", q2$from_cache && file.mtime(q2$qa_file) == m1 && !any(grepl("tmp", list.files(qdir))))
unlink(q1$qa_file); q3 <- water_dist_cached(small, ref, ocean, qdir, qa_plot = TRUE, coastline_id = "ne-test")
check("QA figure: a deleted figure is redrawn from the cached distances (same table)", q3$from_cache && file.exists(q3$qa_file) && isTRUE(all.equal(q3$table$route_km, q1$table$route_km)))
set.seed(1); seed0 <- .Random.seed
invisible(water_dist_cached(small, ref, ocean, tempfile("wdrng"), qa_plot = TRUE, coastline_id = "ne-test"))
check("QA figure: the caller's random-number state is untouched", identical(.Random.seed, seed0))
orig_plot <- water_dist_qa_plot; water_dist_qa_plot <- function(...) stop("boom")
w <- NULL; r <- withCallingHandlers(water_dist_cached(small, ref, ocean, tempfile("wdboom"), qa_plot = TRUE, coastline_id = "ne-test"),
                                    warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
water_dist_qa_plot <- orig_plot
check("QA figure: a plotting failure only warns and the distances are still returned", !is.null(w) && grepl("boom", w) && nrow(r$table) == 3 && is.na(r$qa_file))
cat("Part F - the margin around the study area\n")
kat <- data.frame(lon = 11.5, lat = 57.0)
k <- water_route_km(kat, ref, ocean, res_km = 2.5, pad_km = 50)    # at 2.5 km the northern Oresund is closed: the detour via the Belts needs room
check(sprintf("a margin that cuts the detour is enlarged automatically (%.0f km route; margin 50 -> %g km in %d attempts)",
              k$table$route_km, k$settings$pad_km, k$settings$pad_attempts),
      k$table$status == "ok" && k$settings$pad_attempts > 1 && abs(k$table$route_km - 487) < 6)
cur <- water_route_km(data.frame(lon = 21.1, lat = 55.2), ref, ocean, res_km = 1)   # Curonian Lagoon: enclosed at 1 km
check("an enclosed lagoon is reported unreachable without pointless retries", cur$table$status == "unreachable" && cur$settings$pad_attempts == 1)
cur_p <- water_route_km(data.frame(lon = 21.1, lat = 55.2), ref, ocean, res_km = 1, paths = TRUE)
check("paths = TRUE with no routable site returns NULL instead of crashing", cur_p$table$status == "unreachable" && is.null(cur_p$paths))
k1 <- water_route_km(kat, ref, ocean, res_km = 1, pad_km = 50)
check(sprintf("at 1 km the same site needs no enlargement and gives %.0f km", k1$table$route_km), k1$settings$pad_attempts == 1 && abs(k1$table$route_km - 164) < 4)
cat("Part G - several start points (the nearest entrance wins)\n")
P <- water_dist_parse_points("12.80,55.75; 10.95,55.45")                                    # Oresund; Great Belt
check("water_dist_parse_points reads 'lon,lat; lon,lat' and rejects malformed input",
      identical(dim(P), c(2L, 2L)) && isTRUE(all.equal(P[2, ], c(10.95, 55.45))) && is.null(water_dist_parse_points("")) &&
      inherits(try(water_dist_parse_points("12.8;55.7"), silent = TRUE), "try-error") && inherits(try(water_dist_parse_points("200,55"), silent = TRUE), "try-error"))
g <- data.frame(lon = c(10.5, 13.8, 15.5, 11.5), lat = c(54.5, 54.9, 55.3, 57.0))           # Kiel Bay, Arkona, Bornholm, Kattegat
laea <- "+proj=laea +lat_0=56 +lon_0=13 +datum=WGS84 +units=m +no_defs"                      # ONE fixed grid for all runs below, so they are comparable
m2 <- water_route_km(g, P, ocean, crs = laea, paths = TRUE)
a1 <- water_route_km(g, P[1, ], ocean, crs = laea)$table; a2 <- water_route_km(g, P[2, ], ocean, crs = laea)$table
check("nearest of two start points = the smaller of the two single-start routes", all(abs(m2$table$route_km - pmin(a1$route_km, a2$route_km)) < 0.05))
check("ref_id says which start point was nearer", identical(m2$table$ref_id, ifelse(a1$route_km <= a2$route_km, 1L, 2L)))
check(sprintf("Kiel Bay is reached via the Great Belt: %.0f km (%.0f km from Oresund alone), Arkona stays with Oresund", m2$table$route_km[1], a1$route_km[1]),
      m2$table$ref_id[1] == 2 && m2$table$route_km[1] < 0.6 * a1$route_km[1] && m2$table$ref_id[2] == 1)
check("route polylines run to the nearest start point and have the reported length (+-0.1 %)",
      all(abs(perim(m2$paths) / 1000 / m2$table$route_km[m2$paths$point_row] - 1) < 1e-3))
err <- try(water_route_km(g, rbind(P[1, ], c(14.0, 57.2)), ocean), silent = TRUE)             # 2nd start point is inland Smaland
check("a start point that cannot be placed stops with its number and coordinates",
      inherits(err, "try-error") && grepl("#2", conditionMessage(attr(err, "condition"))) && grepl("14.0000", conditionMessage(attr(err, "condition"))))
gs <- data.frame(site_key = paste(g$lat, g$lon), lon = g$lon, lat = g$lat)
q <- water_dist_cached(gs, P, ocean, tempfile("wdmulti"), negative_side = c(8, 13, 55.95, 60), qa_plot = TRUE, coastline_id = "ne-test")
check("cached helper with two start points: ref_id, signed_route_km (Kattegat negative), figure and a self-describing sidecar",
      all(q$table$ref_id %in% 1:2) && q$table$signed_route_km[4] < 0 && all(q$table$signed_route_km[1:3] > 0) && file.exists(q$qa_file) && file.size(q$qa_file) > 20000 &&
      all(c("route_km", "geodesic_km", "signed_route_km", "ref_id") %in% names(q$meta$columns)) && nrow(q$meta$reference_snapped_lonlat) == 2)

cat("Part I - the HELCOM line as a start box (real coastline)\n")
helcom_box <- water_dist_parse_boxes("lon [10.60, 12.00], lat [57.7395, 57.7415]")                # the entrance to the Baltic Sea (HELCOM): 57 44.43' N
bal <- read.csv(text = "name,lon,lat
Skagerrak,10.90,58.20
Kattegat_N,11.30,57.50
Kattegat,11.50,57.00
Oresund,12.70,55.90
KielBay,10.30,54.60
Bornholm,15.20,55.10
GdanskBay,19.00,54.70
Gotland_E,19.30,57.50
Stockholm,18.80,59.40
Aland,20.00,60.20
GulfFinland,26.00,60.00
RigaBay,23.50,57.80
BothnianSea,19.50,62.50
BothnianBay,23.00,65.00
VistulaLagoon,19.90,54.55
CuronianLagoon,21.10,55.20")
t0 <- Sys.time()
hb <- water_route_km(bal[, c("lon", "lat")], polygons = ocean, ref_box = helcom_box, res_km = 1, paths = TRUE)
cat(sprintf("  HELCOM box, 16 sites, 1 km grid: %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ht <- hb$table; rownames(ht) <- bal$name; print(round(ht[, c("route_km", "geodesic_km", "circuity", "snap_km")], 1)); print(hb$settings$reference_box_info)
bi <- hb$settings$reference_box_info
check(sprintf("the thin HELCOM line (220 m) selects %d passable 1-km cells (of %d it touches; the rest is land: Skagen, the Bohuslan archipelago), spread over %.2f-%.2f E",
              bi$cells_passable, bi$cells_touched, bi$passable_lon_min, bi$passable_lon_max),
      bi$cells_passable >= 60 && bi$cells_touched > bi$cells_passable && bi$passable_lon_min >= 10.6 - 0.02 && bi$passable_lon_max <= 12.0 + 0.02 && bi$passable_lon_max - bi$passable_lon_min > 0.8)
check("Kattegat site 82 km south of the line: route 80-88 km (the line is at 57.7405 N, the site at 57.0 N, 82.4 km), Skagerrak site 51 km north of it: 49-56 km; both positive",
      ht["Kattegat", "route_km"] > 80 && ht["Kattegat", "route_km"] < 88 && ht["Skagerrak", "route_km"] > 49 && ht["Skagerrak", "route_km"] < 56)
check("unsigned: no signed_route_km, all routed distances >= 0, and the same crow-fly distance gives about the same route on both sides of the line",
      !("signed_route_km" %in% names(ht)) && all(ht$route_km[ht$status == "ok"] >= 0) && abs((ht["Skagerrak", "route_km"] / ht["Skagerrak", "geodesic_km"]) - (ht["Kattegat", "route_km"] / ht["Kattegat", "geodesic_km"])) < 0.05)
okh <- ht$status == "ok"
check("every route is at least as long as the crow-fly line to the box (within the cell tolerance), and a Baltic site is far from the entrance",
      all(ht$route_km[okh] >= ht$geodesic_km[okh] - 1.5) && ht["BothnianBay", "route_km"] > 1000 && ht["BothnianBay", "route_km"] > ht["Stockholm", "route_km"])
check("lagoons are flagged exactly as with a start point (Vistula: no passable cell near; Curonian: unreachable)",
      ht["VistulaLagoon", "status"] == "no_passable_cell_near" && ht["CuronianLagoon", "status"] == "unreachable")
check("route polylines have the reported length (+-0.1 %)", all(abs(perim(hb$paths) / 1000 / ht$route_km[hb$paths$point_row] - 1) < 1e-3))
# the box as one start POINT on each of its cells (2 km grid keeps it to ~40 starts): the same numbers by an independent route
sub <- bal[bal$name %in% c("Skagerrak", "Kattegat", "KielBay", "GdanskBay", "Stockholm", "BothnianBay"), ]
laea2 <- "+proj=laea +lat_0=58 +lon_0=15 +datum=WGS84 +units=m +no_defs"
hb2 <- water_route_km(sub[, c("lon", "lat")], polygons = ocean, ref_box = helcom_box, res_km = 2, crs = laea2)
hp2 <- water_route_km(sub[, c("lon", "lat")], ref = hb2$box_cells[[1]], polygons = ocean, res_km = 2, crs = laea2)
check(sprintf("HELCOM box = one start point per cell it stands for (%d cells at 2 km): identical distances (largest difference %.1e km)", nrow(hb2$box_cells[[1]]),
              max(abs(hb2$table$route_km - hp2$table$route_km), na.rm = TRUE)),
      identical(hb2$table$status, hp2$table$status) && max(abs(hb2$table$route_km - hp2$table$route_km), na.rm = TRUE) < 1e-6)
h1 <- ht[match(sub$name, rownames(ht)), "route_km"]; near <- sub$name %in% c("Skagerrak", "Kattegat", "KielBay")
check("1 km and 2 km grids agree within 3 % (or 4 km) on the sites at the entrance (Skagerrak, Kattegat, Kiel Bay)", all(abs(hb2$table$route_km[near] - h1[near]) < pmax(4, 0.03 * h1[near])))
cat(sprintf("  (the Baltic proper is reached through Oresund at 1 km; a 2 km grid closes it, so GdanskBay/Stockholm/BothnianBay come out %s km longer at 2 km: see WATER_DIST_RES_KM)\n",
            paste(round(hb2$table$route_km[!near] - h1[!near]), collapse = "/")))
hb05 <- water_route_km(sub[, c("lon", "lat")], polygons = ocean, ref_box = helcom_box, res_km = 0.5, crs = laea2)
check(sprintf("halving the cell size (1 km -> 0.5 km, %d cells for the line) changes the HELCOM distances by at most 1 %% (or 4 km): largest %.1f km", hb05$settings$reference_box_info$cells_passable,
              max(abs(hb05$table$route_km - h1))), all(abs(hb05$table$route_km - h1) < pmax(4, 0.01 * h1)))
check("EPSG:3035 and the local LAEA agree within 1 % (or 2 km) on the HELCOM distances",
      all(abs(water_route_km(bal[okh, c("lon", "lat")], polygons = ocean, ref_box = helcom_box, res_km = 1, crs = "EPSG:3035")$table$route_km - ht$route_km[okh]) < pmax(2, 0.01 * ht$route_km[okh])))
# optional sign on the same line: Skagerrak (north of it) negative, the magnitudes unchanged
hs <- water_route_km(bal[, c("lon", "lat")], polygons = ocean, ref_box = helcom_box, res_km = 1, negative_side = c(8, 13, 57.75, 60))$table
check("optional sign: the same distances, negated north of the line (Skagerrak negative, Kattegat and the Baltic positive)",
      "signed_route_km" %in% names(hs) && isTRUE(all.equal(hs$route_km, ht$route_km)) && hs$signed_route_km[1] < 0 && all(hs$signed_route_km[hs$status == "ok" & seq_len(nrow(hs)) != 1] > 0))
# cache and figure
bsites <- data.frame(site_key = paste(bal$lat, bal$lon), lon = bal$lon, lat = bal$lat)
hq <- water_dist_cached(bsites, polygons = ocean, cache_dir = tempfile("wdhelcom"), ref_box = helcom_box, qa_plot = TRUE, coastline_id = "ne-test")
check("cache + QA figure for the HELCOM line: figure written, box described in the sidecar, TSV without a sign column",
      file.exists(hq$qa_file) && file.size(hq$qa_file) > 20000 && hq$meta$n_reference_boxes == 1 && hq$meta$reference_box_info$cells_passable == bi$cells_passable &&
      !("signed_route_km" %in% names(hq$table)) && isTRUE(all.equal(hq$table$route_km, ht$route_km)))
invisible(file.copy(hq$qa_file, file.path(Sys.getenv("WD_TEST_FIGURES_DIR", tempdir()), "helcom_box_qa_figure.png"), overwrite = TRUE))
hq2 <- water_dist_cached(bsites, polygons = ocean, cache_dir = tempfile("wdhelcom2"), ref_box = helcom_box, negative_side = c(8, 13, 57.75, 60), qa_plot = TRUE, coastline_id = "ne-test")
check("the QA figure with the optional negative side draws without error", file.exists(hq2$qa_file) && file.size(hq2$qa_file) > 20000)
invisible(file.copy(hq2$qa_file, file.path(Sys.getenv("WD_TEST_FIGURES_DIR", tempdir()), "helcom_box_with_negative_side_qa_figure.png"), overwrite = TRUE))
cat(sprintf("all %d checks passed\n", n_checks))
