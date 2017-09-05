## The function plots the number of dense zones for each threshold value.
## Input: a list of seven elements, as the output of the dense_zones function.
## Output: a plot with the number of dense zones found.

plot_n_zones <- function(zones) {
    if (!is.list(zones))
        stop("'zones' must be an object of type 'list'.")
    plot(zones$zones_count$threshold, zones$zones_count$n_zones, type = "l",
        xlab = "threshold", ylab = "# zones", main =
        paste("Number of dense DNA zones for different values \nof ",
        zones$acctype, " accumulation (threshold)", sep = ""))
    ## finding the point with max slope change.
    d <- diff(diff(zones$zones_count$n_zones))
    max_slope_change_position <- zones$zones_count$threshold[1 +
        which.max(abs(d))]
    points(max_slope_change_position, zones$zones_count$n_zones[which
        (zones$zones_count$threshold == max_slope_change_position)], col =
        "red", lwd = 3)
    legend("topright", legend = c("max slope change"), col = "red", pch = 1)
}
