
#' Create ROV track object
#' 
#' Create a track object from a data.frame. 
#' Preferably, coordinates should be in a projected coordinate system (e.g. UTM).
#' 
#' @name trackROV
#' @details Use \code{create_track} to create a 'track' object and apply filters.
#' Following functions can also be used:
#' \itemize{
#'   \item \code{write_track} - Save coordinates in text file.
#'   \item \code{write_track_st} - Save coordinates in shape file.
#'   \item \code{get_coords} - Get coordinates from a 'track' object.
#'   \item \code{copy_track} - Make a copy of the object.
#'   \item \code{plot} - 3D plot of the track using \code{plotly}.
#'   \item \code{\link{plot_tracks}} - Save tracks of each dive as a html file.
#'   \item \code{\link{set_dives}} - Associate dives based on start and end timestamps.
#'   }
#' Always use \code{copy_track} to make a copy of the track object (e.g. \code{track2 = copy(track1)}),
#' or else filters will be applied to both objects.
#' 
#' Possible filters include:
#' \itemize{
#'   \item \code{\link{filter_depth}} - Filter points based on depth 
#'   \item \code{\link{trim_track}} - Remove points that are far away in the start or end of a track.
#'   \item \code{\link{remove_outliers}} - Remove points that are far away from adjacent points.
#'   \item \code{\link{remove_stopped}} - Detect cloud of points and remove them from the track.
#'   \item \code{\link{distance_filter}} - Distance threshold filter.
#'   \item \code{\link{distance_filter_mean}} - Distance threshold filter followed by mean.
#'   \item \code{\link{rolling_mean}} - Rolling mean filter.
#'   \item \code{\link{simplify_filter}} - Remove points using Douglas-Peuker Algorithm.
#'   \item \code{\link{interpolate_points}} - Interpolate points along track.
#'   \item \code{\link{fill_time_gaps}} - Create a point every second.
#'   }
#' @examples 
#' # Example with dummy track
#' data("simcoords")
#' track_dummy = create_track(simcoords) %>%  
#'     remove_outliers(n = 15, d = 4, select = "Dive1") %>%
#'     rolling_mean(n = 15) %>% 
#'     remove_stopped(n = 61, p = 0.0004, d = 6, t=35)
#' plot(track_dummy)
#' 
#' track_dummy2 = track_dummy %>% copy_track() %>% # use \code{copy_track} or else filters will also affect track_dummy
#'     distance_filter(d = 4)
#' plot(track_dummy2)
#' plot(track_dummy2, select = "Dive2")
#' 
#' write_track(track_dummy2, "track_clean.csv", digits = 2)
#' @param coords A data.frame with the ROV data. See \code{\link{simcoords}} for more details.
#' @param crs A object with coordinate reference system. Check \code{\link[sf]{st_crs}} for options.
#' @param crs_from If supplied, coordinates are converted from \code{crs_from} to \code{crs}.
#' @rdname trackROV
#' @importFrom data.table :=
#' @export
create_track <- function(coords, crs = NA, crs_from = NA) {
    if (!any(class(coords) == "data.frame")) stop("coords must be a data.frame or data.table")
    coords <- data.table::as.data.table(coords)
    
    if (!all(c("Time","Lon","Lat","Depth") %in% colnames(coords))) 
        stop("coords must have 'Time, Lon, Lat, Depth' columns")

    if (!("Dive" %in% colnames(coords))) {
        coords[, Dive := "dive1"]
    }
    
    if (!is.na(crs_from)) {
      x = sf::sf_project(from = sf::st_crs(crs_from), to = sf::st_crs(crs), coords[, c('Lon','Lat', 'Depth')])
      coords[, Lon := x[, 1]]
      coords[, Lat := x[, 2]]
      coords[, Depth := x[, 3]]
    }
    
    data.table::setkey(coords, Dive, Time)
    
    x <- list()
    x$cols <- c("Lon","Lat","Depth")
    x$crs <- sf::st_crs(crs)
    x$interpolated = FALSE
    x$coords_init <- data.table::copy(coords)
    x$coords_init[, removed := "kept"]
    x$coords <- data.table::copy(coords)
    class(x) <- "track"
    return(x)
}


#' @param x 'track' object
#' @param filename Output file name.
#' @param digits  Number of digits to round data
#' @rdname trackROV
#' @export
write_track <- function(x, filename, digits = 2, ...) {
    x$coords[, (x$cols) := lapply(.SD, round, digits = digits), .SDcols=x$cols]
    
    if (length(unique(x$coords$Dive)) == 1) {
      data.table::fwrite(x$coords[, !"Dive"], filename, ...)
    } else {
      data.table::fwrite(x$coords, filename, ...)
    }
}


#' @param as_line logical. If TRUE, save features as lines, else save as points
#' @param ... further arguments passed to \code{\link[data.table]{fwrite}}, \code{\link[sf]{st_write}},
#' \code{\link[plotly]{plot_ly}} or \code{\link[plotly]{plot3d}}.
#' @rdname trackROV
#' @export
write_track_sf <- function(x, filename, as_line = FALSE, ...) {
  tmp <- x$coords %>%
    sf::st_as_sf(coords = x$cols, crs = x$crs)

  if (as_line) { 
    tmp <- tmp %>%
      dplyr::group_by(Dive) %>%
      dplyr::summarise(Start = dplyr::first(Time),
                       End = dplyr::last(Time),
                       geometry = sf::st_combine(geometry)) %>%
      sf::st_cast("LINESTRING")
    
    tmp$Length <- sf::st_length(sf::st_transform(tmp, 4326))
  }

  sf::st_write(tmp, filename, ...)
}


#' @export
print.track <- function(x, ...) {
    cat("track object\n")
    print(x$coords, ...)
}

#' @param data.table logical. TRUE returns a data.table. FALSE returns a data.frame.
#' @rdname trackROV
#' @export
get_coords <- function(x, data.table = FALSE) {
    if (data.table) return(data.table::copy(x$coords)) else return(as.data.frame(x$coords))
}


#' @param select character, indicating which Dive to plot. Default to plot first dive.
#' @param engine character, either "plotly" or "rdl".
#' @param scale_axis logical, scale axis to avoid rounding error in the scene
#' @rdname trackROV
#' @export
plot.track <- function(x, select = NA, engine = "plotly", scale_axis = FALSE, ...) {
    if (is.na(select)) select = x$coords[1, Dive]
    
    dt1 <- x$coords[Dive==select]
    dt2 <- x$coords_init[Dive==select]
    
    if (scale_axis) {
        minimal = min(dt2$Lon)
        dt1[, Lon := (Lon-minimal)*1000]
        dt2[, Lon := (Lon-minimal)*1000]
        
        minimal = min(dt2$Lat)
        dt1[, Lat := (Lat-minimal)*1000]
        dt2[, Lat := (Lat-minimal)*1000]
        
        minimal = max(dt2$Depth)
        dt1[, Depth := (Depth-minimal)*1000]
        dt2[, Depth := (Depth-minimal)*1000]
    }

    if (engine == "plotly") {
        if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' not installed")
        Time_sec = as.numeric(x$coords[Dive==select, Time])
        Time_sec = Time_sec - Time_sec[1]
    
        defaultW <- getOption("warn")
        options(warn = -1) 
        fig <- plotly::plot_ly(type="scatter3d", mode="lines", ...) %>%
            plotly::add_trace(data=dt1, x=~Lon, y=~Lat, z=~Depth, text=~Time,
                      color=Time_sec, line=list(width = 4), name="Track") %>% 
            plotly::add_trace(data=dt2, x = ~Lon, y = ~Lat, z = ~Depth, 
                      color = ~removed, size = 1, mode="markers", hoverinfo="none") %>% 
            plotly::layout(scene = list(aspectmode='data')) %>% 
            plotly::colorbar(title = "Time (sec)")
        
        options(warn = defaultW)
        return(fig)
    }
    else {
        if (!requireNamespace("rgl", quietly = TRUE)) stop("Package 'rgl' not installed")
        cols <- data.table(removed = c("distance_filter", "is_stopped", "kept", "rolling_mean", "outlier"), colors = c(4,6,3,5,2))
        dt2 <- dt2[cols, on = .(removed)]
        rgl::par3d(windowRect = c(100, 100, 1024, 1024))
        rgl::plot3d(dt1$Lon, dt1$Lat, dt1$Depth, type = "l", col = "red", xlab = 'Lon', ylab = 'Lat', zlab = 'Depth', aspect=FALSE, ...)
        rgl::points3d(dt2$Lon, dt2$Lat, dt2$Depth, col = dt2$colors)
        rgl::legend3d("topright", cols$removed, pch = 16, col = cols$colors)
    }
}


#' @rdname trackROV
#' @export
copy_track <- function(x) {
    x$coords_init <- data.table::copy(x$coords_init)
    x$coords <- data.table::copy(x$coords)
    return(x)
}


check_track <- function(x) {
    if(!("track" %in% class(x))) stop("x must be 'track' object. Use create_track first.")
    if (x$interpolated) stop("This filter should not be applied after 'fill_time_gaps'")
}

#' @export
magrittr::`%>%`

utils::globalVariables(c(".", ".SD", "Depth", "Dive", "Lat", "Lon", "Time", "geometry",
                         "i.Depth", "i.Lat", "i.Lon", "id", "keep", "removed"))

#' Sample of ROV track data
#'
#' A data.frame containing a sample of coordinates from two ROV dives.
#'
#' @format A data frame with 4002 rows and 5 variables:
#' \describe{
#'   \item{Time}{timestamp of the track, numeric or POSIXct}
#'   \item{Lon}{Longitude}
#'   \item{Lat}{Latitude}
#'   \item{Depth}{Depth}
#'   \item{Dive}{optional. Character or factor, indicating from which dive is the coordinates}
#' }
"simcoords"
