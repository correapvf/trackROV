# trackROV

This package includes filters to smooth and clear ROV track data,
and possibility GPS data in general. Filters include remove outliers, moving average,
detecte and remove areas where the ROV were stopped, remove points based on a distance thrhold,
and simplify using the Douglas-Peuker Algorithm. The package also provides a plot funtion to
easily visualise the track in a 3D plot.

**Possible filters include:**
* remove_outliers - Remove points that are far away from adjacent points.
* remove_stopped - Detect cloud of points and remove them from the track.
* distance_filter - Distance threshold filter.
* rolling_mean - Rolling mean filter.
* simplify_filter - Remove points using Douglas-Peuker Algorithm.
* interpolate_points - Interpolate points along track.
* fill_time_gaps - Create a point every second.

Check `?trackROV` and `?simcoords` to get started.
## Installation
```r
library(devtools)
install_github("correapvf/trackROV")
```

## Example
```r
data("simcoords")
track_dummy = create_track(simcoords) %>%  
    remove_outliers(n = 15, d = 4, select = "Dive1") %>%
    rolling_mean(n = 15) %>% 
    remove_stopped(n = 61, p = 0.0004, d = 6, t=35)
plot(track_dummy)

track_dummy2 = track_dummy %>% copy_track() %>%
    distance_filter(d = 4)
plot(track_dummy2)
plot(track_dummy2, select = "Dive2")

write_track(track_dummy2, "track_clean.csv", digits = 2)
```