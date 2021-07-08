
# put the time stamp in the n 
join_coords <- function(lon, lat, depth, dive, coords_along) {
    n = length(lon)
    lon.f = rep(0, n)
    lat.f = rep(0, n)
    depth.f = rep(0, n)
    i = 1
    
    lon.n = coords_along[Dive == dive, Lon]
    lat.n = coords_along[Dive == dive, Lat]
    depth.n = coords_along[Dive == dive, Depth]
    n2 = length(lon.n)
    j = 1
    
    while (i <= n) {
        nf = min(j+60, n2)
        dists = sqrt((lon[i] - lon.n[j:nf])^2 + (lat[i] - lat.n[j:nf])^2 + (depth[i] - depth.n[j:nf])^2)
        j = j + which.min(dists) - 1
        
        lon.f[i] = lon.n[j]
        lat.f[i] = lat.n[j]
        depth.f[i] = depth.n[j]
        
        i = i + 1
    }
    
    return(list(lon.f, lat.f, depth.f))
}

#' Simplify points
#' 
#' Remove points using Douglas-Peuker Algorithm.
#' @param x 'track' object.
#' @param dTolerance,preserveTopology Passed to \code{\link[sf]{st_simplify}}.
#' @importFrom magrittr %>%
#' @export
simplify_filter <- function(x, dTolerance = 5, preserveTopology = TRUE){
    check_track(x)
    
    # create a sf object
    tmp <- x$coords %>%
        sf::st_as_sf(coords = x$cols, crs = x$crs) %>%
        dplyr::group_by(Dive) %>%
        dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
        sf::st_cast("LINESTRING")
    
    # simplify
    tmp2 = sf::st_simplify(tmp, preserveTopology=preserveTopology, dTolerance=dTolerance)

    
    # return to data.table
    coords_along <- data.table::as.data.table(sf::st_coordinates(tmp2))
    colnames(coords_along) <- c(x$cols, "Dive")
    coords_along[, Dive := tmp$Dive[Dive]]
    
    x$coords[, keep := FALSE]
    x$coords[coords_along, keep := TRUE, on=.(Dive, Lat, Lon, Depth)]
    
    tmp <- x$coords[keep == FALSE]
    x$coords_init[tmp, removed := "simplify_filter"] # on=.(Dive, Time)

    x$coords = x$coords[keep == TRUE]
    x$coords$keep <- NULL

    return(x)
}
    
    
#' Interpolate points along track
#' 
#' Create points along track and assign original timestamp to the new points.
#' @param x 'track' object.
#' @param density Number of points per distance unit. Passed to \code{\link[sf]{st_line_sample}}.
#' @export
interpolate_points <- function(x, density = 2){
    check_track(x)
    
    # create a sf object
    tmp <- x$coords %>%
        sf::st_as_sf(coords = x$cols, crs = x$crs) %>%
        dplyr::group_by(Dive) %>%
        dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
        sf::st_cast("LINESTRING")
    
    # create point along line
    tmp2 = sf::st_line_sample(tmp, density=density)
    
    # return to data.table
    coords_along <- data.table::as.data.table(sf::st_coordinates(tmp2))
    colnames(coords_along) <- c(x$cols, "Dive")
    coords_along[, Dive := tmp$Dive[Dive]]
    
    # set original points to their new locations
    x$coords = x$coords_init[removed %in% c("simplify_filter","distance_filter","is_stopped","kept")]
    x$coords[, (x$cols) := lapply(.SD, data.table::frollmean, n=7, align = "center"), by=.(Dive), .SDcols=x$cols]
    x$coords = na.omit(x$coords)
    
    x$coords[, (x$cols) := join_coords(Lon, Lat, Depth, data.table::first(Dive), coords_along), by=.(Dive)]
    
    return(x)
}


#' Create a point every second
#' 
#' Fill gaps in the timestamp of track, generating a point every second.
#' Missing points are interpolated based on \code{type}.
#' @param x 'track' object.
#' @param type interpolation to be used. One of 'locf', 'approx' or 'spline'.
#' @note This should be the last filter used in the track data.
#' @export
fill_time_gaps <- function(x, type = "locf") {
    check_track(x)
    
    # interpolate points for each second
    coords = x$coords[, .(Time = seq(data.table::first(Time), data.table::last(Time), by="sec")), by=.(Dive)]
    coords = x$coords[coords,.(Dive, Time, Lon, Lat, Depth), on=.(Dive, Time)]
    
    if (type == "locf") {
        data.table::setnafill(coords, type="locf", cols=x$cols)
    } else if (type == "approx") {
        coords[, (x$cols) := lapply(.SD, zoo::na.approx), by=.(Dive), .SDcols=x$cols]
    } else if (type == "spline") {
        coords[, (x$cols) := lapply(.SD, zoo::na.spline), by=.(Dive), .SDcols=x$cols]
    } else {
        stop("type should be one of 'locf', 'approx' or 'spline'")
    }
    
    x$coords <- coords
    x$interpolated = TRUE
    return(x)
}
