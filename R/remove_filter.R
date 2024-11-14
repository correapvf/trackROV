
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
    
    n_tail <- ceiling(n/3)
    n_center <- ceiling(n/2)
    tmp <- x$coords[, lapply(.SD, data.table::frollapply, n=n, align = "center", fill=0,
                            FUN= function(x) x[n_center] - stats::median(x[c(1:n_tail, (n-n_tail+1):n)])),
                   by = .(Dive), .SDcols = x$cols]
    
    x$coords$Lon[abs(tmp$Lon) >= d] <- NA
    x$coords$Lat[abs(tmp$Lat) >= d] <- NA
    x$coords$Depth[abs(tmp$Depth) >= d] <- NA
    
    tmp <- na.omit(x$coords, x$cols, invert=TRUE)
    x$coords_init[tmp, removed := "outlier"] # on=.(Dive, Time)
    
    if (!is.null(select)) {
        coords[x$coords, (x$cols) := .(i.Lon, i.Lat, i.Depth)]
        x$coords <- coords
    }
    
    x$coords <- na.omit(x$coords, x$cols)
    
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
#' @param include_mean bool. Just remove the points detected as stopped (default), or calculate the mean and add in the track.
#' @export
remove_stopped <- function(x, n = 61, p = 0.0001, d = 6, t=30, include_mean=FALSE, select = NULL){
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
    
    tmpc <- list()
    
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
            
            tmpmean <- x$coords[tmpi, lapply(.SD, mean), .SDcols = -c('Dive')]
            tmpmean$Dive <- x$coords[tmpi[1], Dive]
            tmpc <- rbind(tmpc, tmpmean)
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
        x$coords <- x$coords[tmp$keep == TRUE]
    }
    
    if (include_mean) {
      x$coords <- rbind(x$coords, tmpc)
      tmpc$removed <- "stopped_mean"
      x$coords_init <- rbind(x$coords_init, tmpc)
      data.table::setkey(x$coords, Dive, Time)
      data.table::setkey(x$coords_init, Dive, Time)
    }
    
    return(x)
}


#' Filter points based on depth 
#' 
#' Remove points where depth is greater than \code{depth}.
#' @param x 'track' object.
#' @param depth threshold to filter points, same units original data.
#' Points with depth above this value will be removed.
#' @param select character. The filter is applied only in these dives.
#' @export
filter_depth <- function(x, depth = -5, select = NULL) {
  check_track(x)
  
  if (!is.null(select)) {
    x$coords <- x$coords[(Depth <= depth) & (Dive %in% select)]
    x$coords_init[(Depth > depth) & (Dive %in% select), removed := "shallow"]
  } else {
    x$coords <- x$coords[Depth <= depth]
    x$coords_init[Depth > depth, removed := "shallow"]
  }
  
  return(x)
}


#' Trim start and end of each track 
#' 
#' Remove points that are far away in the start or end of a track.
#' 
#' This filter will analyze only the ends of each track. If distance is greater than \code{d},
#' any previous/following points will be removed.
#' @param x 'track' object.
#' @param d threshold to remove outliers, same units original data.
#' @param nend maximum number of points at each end that can removed.
#' @param select character. The filter is applied only in these dives.
#' @export
trim_track <- function(x, d = 10, nend = 100, select = NULL) {
  check_track(x)
  
  if (!is.null(select)) {
    x$coords[Dive %in% select, keep := trim_helper(Lat, Lon, Depth, d=d, nend=nend), by = .(Dive)]
  } else {
    x$coords[, keep := trim_helper(Lat, Lon, Depth, d=d, nend=nend), by = .(Dive)]
  }
  
  tmp2 <- x$coords[keep == FALSE]
  x$coords_init[tmp2, removed := "trimmed"] # on=.(Dive, Time)
  
  x$coords <- x$coords[keep == TRUE]
  x$coords$keep <- NULL
  
}


# Function used for trim_track
trim_helper <- function(Lat, Lon, Depth, d, nend) {
  track <- data.table(Lat=Lat, Lon=Lon, Depth=Depth)
  tmp <- head(track, nend)
  tmp[, dist := sqrt((Lat - shift(Lat, -1))^2 + (Lon - shift(Lon, -1))^2 + (Depth - shift(Depth, -1))^2)]
  istart <- tail(which(tmp$dist > d), 1) + 1
  if (length(istart) == 0) istart <- 1
  
  tmp <- tail(track, nend)
  tmp[, dist := sqrt((Lat - shift(Lat))^2 + (Lon - shift(Lon))^2 + (Depth - shift(Depth))^2)]
  iend <- nrow(track) - (nend - head(which(tmp$dist > d), 1)) - 1
  if (length(iend) == 0) iend <- nrow(track)
  
  n <- nend - 1
  if (isart == n | iend == n) warning("'trim_track' removed all points from one end.
                                 Consider reviewing and raising 'd' or 'nend'")
  
  track[, keep := FALSE]
  track[istart:iend, keep := TRUE]
  return(track$keep)
}
