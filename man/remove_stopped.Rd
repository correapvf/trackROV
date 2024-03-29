% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_filter.R
\name{remove_stopped}
\alias{remove_stopped}
\title{Remove points were ROV is stopped}
\usage{
remove_stopped(x, n = 61, p = 1e-04, d = 6, t = 30, select = NULL)
}
\arguments{
\item{x}{'track' object.}

\item{n}{window size to calculate proximity between points.}

\item{p}{threshold to detect were ROV is stopped. Lower values, more likely to detect stops.}

\item{d}{distance buffer around detected points, same units original data. All points within this buffer are removed.}

\item{t}{'time' buffer around detected points, in seconds. All points that within this period are removed.}

\item{select}{character. The filter is applied only in these dives.}
}
\description{
Detect cloud of points and remove them from the track.
}
\details{
This filter will calculate the inverse distance from \code{n} adjacent points.
Than, it use \code{p} to detect cloud of points, with are indicative of places were the ROV is stopped
and there is too much noise in the track. Setting \code{p} will require some trial-and-error.
Than, it creates a buffer (both in distance and time) for each cloud of points that was detected, 
and all points within this buffer are removed.
}
