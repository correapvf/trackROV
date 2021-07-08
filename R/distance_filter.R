
# Distance threshold filter
dist_thr <- function(lon, lat, depth, d) {
    result = rep(FALSE, length(lon))
    i = 1
    final = length(lon)
    result[1] = TRUE
    j = 2
    
    while (i <= final) {
        dists = sqrt((lon[i] - lon[j])^2 + (lat[i] - lat[j])^2 + (depth[i] - depth[j])^2)
        
        if (dists >= d) {
            result[i] = TRUE
            j = i
        }
        
        i = i + 1
    }
    
    return(result)
}

#' Distance threshold filter
#' 
#' Only keep points after a distance \code{d} has been reached.
#' Check \href{https://www.gpsvisualizer.com/tutorials/track_filters.html}{here} for a more detailed
#' explanation of how the filter works.
#' @param x 'track' object.
#' @param d Distance threshold, same units original data.
#' @param select character. The filter is applied only in these dives.
#' @export
distance_filter <- function(x, d = 5, select = NULL) {
    check_track(x)
    
    if (is.null(select)) {
        x$coords[, keep := dist_thr(Lon, Lat, Depth, d), by=.(Dive)]
    } else {
        x$coords[, keep := TRUE]
        x$coords[Dive %in% select, keep := dist_thr(Lon, Lat, Depth, d), by=.(Dive)]
    }
    
    tmp <- x$coords[keep == FALSE]
    x$coords_init[tmp, removed := "distance_filter"] # on=.(Dive, Time)
    
    x$coords = x$coords[keep == TRUE]
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
rolling_mean <- function(x, n = 15, select = NULL){
    if(!("track" %in% class(x))) stop("x must be 'track' object. Use create_track first.")
    
    if (is.null(select)) {
        x$coords[, (x$cols) := lapply(.SD, data.table::frollmean, n=n, align = "center"), by=.(Dive), .SDcols=x$cols]
    } else {
        x$coords[Dive %in% select, (x$cols) := lapply(.SD, data.table::frollmean, n=n, align = "center"), by=.(Dive), .SDcols=x$cols]
    }
    
    tmp <- na.omit(x$coords, x$cols, invert=TRUE)
    x$coords_init[tmp, removed := "rolling_mean"] # on=.(Dive, Time)
    
    x$coords = na.omit(x$coords, x$cols)
    
    return(x)
}
