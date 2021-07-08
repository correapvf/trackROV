
#' Remove outliers
#' 
#' Remove points that are far away from adjacent points.
#' 
#' This filter will divide the window in three parts. A median is calculated from both ends,
#' and if the center value is more distant than \code{d} from the mean, this value is considered an outlier.
#' This process is repeated for each coordinate.
#' @param x 'track' object.
#' @param n window size to search for outlier.
#' @param d threshold to remove outliers, same units original data.
#' @param select character. The filter is applied only in these dives.
#' @importFrom stats na.omit
#' @export
remove_outliers <- function(x, n = 7, d = 5, select = NULL) {
    check_track(x)
    
    if (n < 3) stop("n must be greater or equal to 3")
    
    if (!is.null(select)) {
        coords <- data.table::copy(x$coords)
        x$coords <- x$coords[Dive %in% select]
    }
    
    n_tail = ceiling(n/3)
    n_center = ceiling(n/2)
    tmp = x$coords[, lapply(.SD, data.table::frollapply, n=n, align = "center", fill=0,
                            FUN= function(x) x[n_center] - stats::median(x[c(1:n_tail, (n-n_tail+1):n)])),
                   by=.(Dive), .SDcols=x$cols]
    
    x$coords$Lon[abs(tmp$Lon) >= d] <- NA
    x$coords$Lat[abs(tmp$Lat) >= d] <- NA
    x$coords$Depth[abs(tmp$Depth) >= d] <- NA
    
    tmp <- na.omit(x$coords, x$cols, invert=TRUE)
    x$coords_init[tmp, removed := "outlier"] # on=.(Dive, Time)
    
    if (!is.null(select)) {
        coords[x$coords, (x$cols) := .(i.Lon, i.Lat, i.Depth)]
        x$coords <- coords
    }
    
    x$coords = na.omit(x$coords, x$cols)
    
    return(x)
}



#' Remove points were ROV is stopped
#'
#' Detect cloud of points and remove them from the track.
#' 
#' This filter will calculate the inverse distance from \code{n} adjacent points.
#' Than, it use \code{p} to detect cloud of points, with are indicative of places were the ROV is stopped
#' and there is too much noise in the track. Setting \code{p} will require some trial-and-error.
#' Than, it creates a buffer (both in distance and time) for each cloud of points that was detected, 
#' and all points within this buffer are removed.
#' @param x 'track' object.
#' @param n window size to calculate proximity between points.
#' @param p threshold to detect were ROV is stopped. Lower values, more likely to detect stops.
#' @param d distance buffer around detected points, same units original data. All points within this buffer are removed.
#' @param t 'time' buffer around detected points, in seconds. All points that within this period are removed.
#' @param select character. The filter is applied only in these dives.
#' @export
remove_stopped <- function(x, n = 61, p = 0.0001, d = 6, t=30, select = NULL){
    check_track(x)
    
    if (!is.null(select)) {
        coords <- data.table::copy(x$coords)
        x$coords <- x$coords[Dive %in% select]
    }
    
    density = x$coords[, lapply(.SD, data.table::frollapply, n=n, align = "center", fill=0,
                                FUN= function(x) 1/sum(stats::dist(x)) ), by=.(Dive), .SDcols=x$cols]
    
    tmp = density[, lapply(.SD, `>`, y=p), .SDcols=x$cols]
    tmp[, all := apply(tmp, 1, all)]
    tmp[, id := data.table::rleid(all)]
    tmp[, keep := TRUE]
    
    for (i in 1:max(tmp$id)) {
        index = which(tmp$id == i)
        if (tmp$all[index[1]]) {
            xlim = range(x$coords[index, Lon]) + c(-d, d)
            ylim = range(x$coords[index, Lat]) + c(-d, d)
            zlim = range(x$coords[index, Depth]) + c(-d, d)
            mlim = range(x$coords[index, Time]) + c(-t, t)
            
            tmpi = with(x$coords, which((Lon > xlim[1] & Lon < xlim[2] &
                                         Lat > ylim[1] & Lat < ylim[2] &
                                         Depth > zlim[1] & Depth < zlim[2]) |
                                         (Time > mlim[1] & Time < mlim[2])
                                        ))
            
            tmp[tmpi, keep := FALSE]
        }
    }
    
    tmp2 <- x$coords[tmp$keep == FALSE]
    x$coords_init[tmp2, removed := "is_stopped"] # on=.(Dive, Time)
    
    if (!is.null(select)) {
        coords[, keep := TRUE]
        coords[Dive %in% select, keep := tmp$keep]
        x$coords <- coords[keep == TRUE]
        x$coords$keep <- NULL
    } else {
        x$coords = x$coords[tmp$keep == TRUE]
    }
    
    return(x)
}
