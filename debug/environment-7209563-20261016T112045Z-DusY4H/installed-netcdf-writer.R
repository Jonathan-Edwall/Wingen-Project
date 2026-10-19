function (x, file, missval = -999, prec = "double", ...) 
{
    if (!inherits(x, "stars")) 
        stop("Expect `stars` but got", class(x))
    if (requireNamespace("ncdf4") && requireNamespace("units")) {
        dms <- stars::st_dimensions(x)
        dim_nv <- ncdf4::ncdim_def(name = "nv", units = "", vals = 1:2, 
            create_dimvar = FALSE)
        nc_dims <- list()
        nc_vars <- list()
        dim_vals <- list()
        for (dm in dimnames(x)) {
            dim_vals[[dm]] <- unname(as.matrix(suppressMessages(dplyr::bind_cols(lapply(c("center", 
                "start", "end"), function(w) stars::st_get_dimension_values(x, 
                dm, where = w))))))
            unit <- dms[[dm]]$refsys
            unit <- ifelse(is.character(unit), unit, NA_character_)
            if (dm == "elevation") 
                unit <- "m"
            if (dm == "longitude") 
                unit <- "degrees_east"
            if (dm == "latitude") 
                unit <- "degrees_north"
            if (dm == "time") {
                dim_vals[[dm]] <- matrix(as.numeric(lubridate::seconds(lubridate::as_datetime(dim_vals[[dm]]))), 
                  nrow = nrow(dim_vals[[dm]]))
                unit <- "seconds since 1970-01-01 00:00:00 UTC"
            }
            nc_dims[[dm]] <- ncdf4::ncdim_def(name = dm, units = unit, 
                vals = dim_vals[[dm]][, 1])
            nc_vars[[sprintf("var_%s_bnds", dm)]] <- ncdf4::ncvar_def(paste(dm, 
                "bnds", sep = "_"), units = unit, dim = list(dim_nv, 
                nc_dims[[dm]]), prec = prec)
        }
        target_dims <- c("longitude", "latitude", "elevation", 
            "time")
        matched_dims <- target_dims[target_dims %in% dimnames(x)]
        for (v_nm in names(x)) {
            nc_vars[[v_nm]] <- ncdf4::ncvar_def(name = v_nm, 
                units = units::deparse_unit(x[[v_nm]]), dim = nc_dims[matched_dims], 
                missval = missval, prec = prec)
        }
        nc_vars[["crs_var"]] <- ncdf4::ncvar_def(name = "crs", 
            units = "", dim = list(), prec = "char")
        nc_file <- ncdf4::nc_create(file, nc_vars)
        for (dm in dimnames(x)) {
            ncdf4::ncatt_put(nc_file, varid = dm, attname = "bounds", 
                attval = sprintf("var_%s_bnds", dm))
            ncdf4::ncvar_put(nc_file, varid = nc_vars[[sprintf("var_%s_bnds", 
                dm)]], vals = t(dim_vals[[dm]][, -1]))
        }
        ncdf4::ncatt_put(nc_file, "crs", "grid_mapping_name", 
            "latitude_longitude")
        ncdf4::ncatt_put(nc_file, "crs", "crs_wkt", sf::st_crs(x)$wkt)
        for (v_nm in names(x)) {
            ncdf4::ncvar_put(nc_file, varid = nc_vars[[v_nm]], 
                vals = as.array(aperm(x[v_nm], matched_dims)[[v_nm]]))
            ncdf4::ncatt_put(nc_file, v_nm, "grid_mapping", "crs")
        }
        ncdf4::nc_close(nc_file)
    }
}
