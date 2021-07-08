
#' Create ROV track object
#' 
#' Create a track object from a data.frame. 
#' Preferably, coordinates should be in a projected coordinate system (e.g. UTM).
#' 
#' @name trackROV
#' @details Use \code{create_track} to create a 'track' object and apply filters.
#' Following functions can also be used:
#' \itemize{
#'   \item \code{write_track} - Save coordinates in file.
#'   \item \code{get_coords} - Get coordinates from a 'track' object.
#'   \item \code{copy_track} - Make a copy of the object.
#'   \item \code{plot} - 3D plot of the track using \code{plotly}.
#'   }
#' Always use \code{copy_track} to make a copy of the track object (e.g. \code{track2 = copy(track1)}),
#' or else filters will be applied to both objects.
#' 
#' Possible filters include:
#' \itemize{
#'   \item \code{\link{remove_outliers}} - Remove points that are far away from adjacent points.
#'   \item \code{\link{remove_stopped}} - Detect cloud of points and remove them from the track.
#'   \item \code{\link{distance_filter}} - Distance threshold filter.
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
#' track_dummy2 = track_dummy2 %>% copy_track() %>% # use \code{copy_track} or else will also affect track_dummy
#'     distance_filter(d = 4)
#' plot(track_dummy2)
#' plot(track_dummy2, select = "Dive2")
#' 
#' write_track(track_dummy2, "track_clean.csv", digits = 2)
#' @param coords A data.frame with the ROV data. See \code{\link{simcoords}} for more details.
#' @param crs A object with coordinate reference system. Check \code{\link[sf]{st_crs}} for options.
#' @rdname trackROV
#' @importFrom data.table :=
#' @export
create_track <- function(coords, crs = NA) {
    if (!any(class(coords) == "data.frame")) stop("coords must be a data.frame or data.table")
    data.table::setDT(coords)
    
    if (!all(colnames(coords) %in% c("Dive","Time","Lon","Lat","Depth"))) 
        stop("coords must have '[Dive,] Time, Lon, Lat, Depth' columns")

    if (!("Dive" %in% colnames(coords))) {
        coords[, Dive := "dive"]
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
#' @param ... further arguments passed to \code{\link[data.table]{fwrite}} or \code{\link[plotly]{plot_ly}}.
#' @rdname trackROV
#' @export
write_track <- function(x, filename, digits = 2, ...) {
    x$coords[, (x$cols) := lapply(.SD, round, digits = digits), .SDcols=x$cols]
    data.table::fwrite(x$coords, filename, ...)
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
#' @rdname trackROV
#' @export
plot.track <- function(x, select = NA, ...) {
    if (is.na(select)) select = x$coords[1, Dive]
    Time_sec = as.numeric(x$coords[Dive==select, Time])
    Time_sec = Time_sec - Time_sec[1]
    
    defaultW <- getOption("warn")
    options(warn = -1) 
    fig <- plotly::plot_ly(type="scatter3d", mode="lines", ...) %>%
        plotly::add_trace(data=x$coords[Dive==select], x=~Lon, y=~Lat, z=~Depth, text=~Time,
                  color=Time_sec, line=list(width = 4), name="Track") %>% 
        plotly::add_trace(data=x$coords_init[Dive==select], x = ~Lon, y = ~Lat, z = ~Depth, 
                  color = ~removed, size = 1, mode="markers", hoverinfo="none") %>% 
        plotly::layout(scene = list(aspectmode='data')) %>% 
        plotly::colorbar(title = "Time (sec)")
    
    options(warn = defaultW)
    return(fig)
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
    if (x$interpolated) stop("This filter should not be applied after 'fill_by_time'")
}

#' @export
magrittr::`%>%`

utils::globalVariables(c("." ,".SD", "Depth", "Dive", "Lat", "Lon", "Time", "geometry",
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
