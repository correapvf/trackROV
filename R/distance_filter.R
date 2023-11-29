
# Distance threshold filter
dist_thr <- function(lon, lat, depth, d) {
    final <- length(lon)
    result <- rep(FALSE, final)
    i <- 1L
    result[1L] <- TRUE
    j <- 2L

    while (i <= final) {
        dists <- sqrt((lon[i] - lon[j])^2 + (lat[i] - lat[j])^2 + (depth[i] - depth[j])^2)

        if (dists >= d) {
            result[i] <- TRUE
            j <- i
        }

        i <- i + 1L
    }

    return(result)
}

#' Distance threshold filter
#'
#' Only keep points after a distance \code{d} has been reached.
#' Check \href{https://www.gpsvisualizer.com/tutorials/track_filters.html}{here} for a more detailed
#' explanation of how the filter works.
#' 
#' \code{distance_filter_mean} is similar, but it adds an extra step where the mean point is calculated
#' for each group based on the distance threshold. This can minimize the effect of outliers.
#' @param x 'track' object.
#' @param d Distance threshold, same units original data.
#' @param select character. The filter is applied only in these dives.
#' @aliases distance_filter_mean
#' @export
distance_filter <- function(x, d = 5, select = NULL) {
    check_track(x)

    if (is.null(select)) {
        x$coords[, keep := dist_thr(Lon, Lat, Depth, d), by = Dive]
    } else {
        x$coords[, keep := TRUE]
        x$coords[Dive %in% select, keep := dist_thr(Lon, Lat, Depth, d), by = Dive]
    }

    x$coords[, keep := replace(keep, .N, TRUE), by = Dive] # keep last coordinate

    tmp <- x$coords[keep == FALSE]
    x$coords_init[tmp, removed := "distance_filter"] # on=.(Dive, Time)

    x$coords <- x$coords[keep == TRUE]
    x$coords$keep <- NULL

    return(x)
}



#' Rolling mean filter
#'
#' Apply a simple rolling mean filter.
#' @param x 'track' object.
#' @param n window size to calculate the mean.
#' @param select character. The filter is applied only in these dives.
#' @export
rolling_mean <- function(x, n = 15, select = NULL) {
    if(!("track" %in% class(x))) stop("x must be 'track' object. Use create_track first.")

    if (is.null(select)) {
        x$coords[, (x$cols) := lapply(.SD, data.table::frollmean, n = n, align = "center"), by = .(Dive), .SDcols = x$cols]
    } else {
        x$coords[Dive %in% select, (x$cols) := lapply(.SD, data.table::frollmean, n=n, align = "center"), by=.(Dive), .SDcols=x$cols]
    }

    tmp <- na.omit(x$coords, x$cols, invert = TRUE)
    x$coords_init[tmp, removed := "rolling_mean"] # on=.(Dive, Time)

    x$coords <- na.omit(x$coords, x$cols)

    return(x)
}


#' @rdname distance_filter
#' @export
distance_filter_mean <- function(x, d = 5, select = NULL) {
    check_track(x)

    if (is.null(select)) {
        x$coords[, group := dist_thr_mean(Lon, Lat, Depth, d), by=.(Dive)]
    } else {
        x$coords[, group := 0L]
        x$coords[Dive %in% select, keep := dist_thr_mean(Lon, Lat, Depth, d), by = .(Dive)]
    }

    result <- x$coords[, c(lapply(.SD, mean), .(Time = quantile(Time, probs = 0.5, type = 3))),
        by = .(Dive, group), .SDcols = x$cols]

    x$coords[, keep := FALSE]
    x$coords[result, keep := TRUE, on = .(Dive, Time)]

    tmp <- x$coords[keep == FALSE]
    x$coords_init[tmp, removed := "dist_filter_mean"]

      x$coords <- result
    x$coords$group <- NULL

  return(x)
}


# Distance threshold filter
dist_thr_mean <- function(lon, lat, depth, d) {
    final <- length(lon)
    result <- rep(0L, final)
    i <- 2L
    j <- 1L
    dists <- 0
    group <- 1L
    result[1L] <- 1L

    # first point
    while (dists <= d / 2) {
        dists <- sqrt((lon[i] - lon[j])^2 + (lat[i] - lat[j])^2 + (depth[i] - depth[j])^2)
        result[i] <- group
        i <- i + 1L
    }

    group <- group + 1L
    j <- i - 1L

    while (i <= final) {
        dists <- sqrt((lon[i] - lon[j])^2 + (lat[i] - lat[j])^2 + (depth[i] - depth[j])^2)

        if (dists >= d) {
            group <- group + 1L
            j <- i
        }
        result[i] <- group
        i <- i + 1L
    }

    return(result)
    }