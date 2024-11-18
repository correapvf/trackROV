
#' Save plots of all dives
#' 
#' Plot and save tracks of each dive in a track object as a html file.
#' 
#' Only the \code{plotly} engine is supported.
#' @param x 'track' object.
#' @param ... further arguments passed to \code{plot.track}.
#' @export
save_tracks <- function(x, ...) {
  dives <- unique(x$coords$Dive)
  for (dive in dives) {
    fig <- plot.track(x, select = dive, ...)
    htmlwidgets::saveWidget(plotly::as_widget(fig), paste0(dive, ".html"))
    print("Saved file ", dive, ".html!")
  }
}


#' Generate Dive column
#' 
#' Associate dives based on start and end timestamps.
#' 
#' This function will overwrite the dives column based on timestamps of start and end.
#' Any point not contained in the interval of the new dives will
#' be removed. This function is useful when working with tracks that are not separated by
#' dives originally, of if you wish to crop a dive in multiple transects.
#' @param x 'track' object.
#' @param dive_ranges A data.frame or data.table with columns:
#' \itemize{
#'   \item \code{Dive} - Name of the dive.
#'   \item \code{start} - timestamp of the start of the dive, same format as the \code{Time}
#'    used to create the track.
#'   \item \code{end} - timestamp of the end of the dive, same format as \code{start}.
#'   }
#' @export
set_dives <- function(x, dive_ranges) {
  time_ranges <- data.table::as.data.table(dive_ranges)
  data.table::setkey(time_ranges, start, end)
  x$coords[, dummy := Time]
  x$coords <- data.table::foverlaps(x$coords, time_ranges, by.x=c("Time","dummy"))
  x$coords[, c("dummy", "i.Dive", "start", "end") := NULL]
  x$coords <- x$coords[!is.na(Dive)]
  data.table::setkey(x$coords, Dive, Time)
  
  x$coords_init <- data.table::copy(x$coords)
  x$coords_init[, removed := "kept"]
  return(x)
}
