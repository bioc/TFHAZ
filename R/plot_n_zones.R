## The function plots the number of dense zones for each threshold value.
## Input:
##    zones: a list of eight elements, as the output of the dense_zones
##    function.
##    chr: optional argument, needed if the input is found with chr = "all";
##    a string or a vector containing strings with the name of the chromosome(s)
##    (e.g., "chr1" or c("chr1", "chr4"))
## Output: a plot(or more) with the number of dense zones found.

plot_n_zones <- function(zones, chr = NULL) {
    if (!is.list(zones))
        stop("'zones' must be an object of type 'list'.")
    if (zones$chr != "all") {
        chr <- zones$chr
        }
    for (i in seq_along(chr)) {
    plot(eval(parse(text = paste("zones$zones_count$", chr[i],"$threshold",
                                sep = ""))),
        eval(parse(text = paste("zones$zones_count$", chr[i],"$n_zones",
                                sep = ""))), type = "l",
        xlab = "threshold", ylab = "# zones", main =
        paste(chr[i],": number of dense DNA zones for different values \nof ",
        zones$acctype, " accumulation (threshold)", sep = ""))
    ## finding the point with max slope change.
    d <- diff(diff(eval(parse(text = paste("zones$zones_count$", chr[i],
                                        "$n_zones", sep = "")))))
    max_slope_change_position <- eval(parse(text = paste("zones$zones_count$",
                                chr[i],"$threshold", sep = "")))[1 +
                                which.max(abs(d))]
    points(max_slope_change_position, eval(parse(text =
                                                paste("zones$zones_count$",
                                                chr[i],"$n_zones",
                                                sep = "")))[which
        (eval(parse(text = paste("zones$zones_count$", chr[i],"$threshold",
                            sep = ""))) == max_slope_change_position)],
                            col = "red", lwd = 3)
    legend("topright", legend = c("max slope change"), col = "red", pch = 1)

    }
}
