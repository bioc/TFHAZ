## This function finds transcription factor high accumulation DNA zones (TFHAZ).
## Each TFHAZ is formed by contiguous bases with accumulation higher than the
## threshold (TH).
## input:
##    accumulation: list of four elements containing: a sparse vector with
##    accumulation values (e.g.,obtained with the accumulation function),
##    the accumulation type, a chromosome name, and the half-width of the
##    window used for the accumulation count.
## output: A list of eight elements:
##    n_zones: an integer containing the number of high accumulation zones
##    obtained.
##    n_bases: an integer containing the total number of bases belonging to
##    the accumulation zones obtained.
##    lengths: a vector containing the considered threshold value and min,
##    max, mean, median and standard deviation of the accumulation
##    zone lengths obtained.
##    distances: a vector containing the considered threshold value and  min,
##    max, mean, median and standard deviation of the distances between adjacent
##    accumulation zones obtained.
##    TH: a number with the threshold value found.
##    acctype: a string with the accumulation type used.
##    chr: a string with the chromosome name associated with the
##    accumulation vector used.
##    w: an integer with half-width of the window used to calculate the
##    accumulation vector.
## A ".bed" file with the chromosome and genomic coordinates of the accumulation
## zones found is created.
## A ".png" file with the plot of the TFHAZ found along the cromosome is
## created.


high_accumulation_zones <- function(accumulation) {

    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")
    if (class(accumulation$accvector) != "Rle")
        stop("'accvector' element of the list 'accumulation' must be of class
        'Rle'.")

    ## input accumulation vector
    acc_tf <- accumulation$accvector

    ## finding the threshold
    TH <- mean(acc_tf[acc_tf > 0]) + 2 * sd(acc_tf[acc_tf > 0])
    ## finding high accumulation zones
    zones <- slice(acc_tf, TH)
    ## finding the number of bases belonging to the zones
    n_bases <- sum(width(zones))
    ## finding the number of zones
    n_zones <- length(zones)
    ## finding the bases belonging to the zones
    bases <- Rle()
    for (i in seq_along(zones)) {
        bases <- append(bases, start(zones)[i]:end(zones)[i])
    }

    ## plot
    png(filename = paste("high_accumulation_zones_TH_", round(TH, digits = 1),
        "_", accumulation$acctype, "_acc_w_", accumulation$w, "_",
        accumulation$chr, ".png", sep = ""))
    plot(accumulation$accvector, type = "l", xlab = "base", ylab = paste
        ("# of ", accumulation$acctype, "s", sep = ""))
    points(bases, rep(0, length(bases)), col = "red", pch = 15)
    abline(h = TH, col = "red")
    legend("topleft", legend = c("threshold", "high accumulation zones"),
        col = "red", pch = c(NA, 15), lty = c(1, NA))
    dev.off()

    ## writing on files chromosome and positions of starting and ending points
    write.table(data.frame(rep(accumulation$chr, length(zones)), start(zones)
        - 1, end(zones)), file = paste(accumulation$acctype, "_acc_w_",
        accumulation$w, "_", accumulation$chr, "_dense_zones_th_",
        round(TH, digits = 1), ".bed", sep = ""), row.names = FALSE,
        col.names = FALSE, quote = FALSE, sep = "\t")

    ## finding elements of length dataframe
    length_zone_min <- min(width(zones))
    length_zone_max <- max(width(zones))
    length_zone_mean <- mean(width(zones))
    length_zone_median <- median(width(zones))
    length_zone_sd <- sd(width(zones))

    ## calculating distances of dense zones
    if (length(zones) > 1) {

        dist_zone <- start(zones)[-1] - end(zones[-length(end(zones))])
        dist_zone_min <- min(dist_zone)
        dist_zone_max <- max(dist_zone)
        dist_zone_mean <- mean(dist_zone)
        dist_zone_median <- median(dist_zone)
        dist_zone_sd <- sd(dist_zone)
    } else {
        dist_zone <- NA
        dist_zone_min <- NA
        dist_zone_max <- NA
        dist_zone_mean <- NA
        dist_zone_median <- NA
        dist_zone_sd <- NA
    }

    ## vectors creation
    lengths <- c(TH, n_zones, length_zone_min, length_zone_max,
        length_zone_mean, length_zone_median, length_zone_sd)
    names(lengths) = c("TH", "n_zones", "length_zone_min", "length_zone_max",
        "length_zone_mean", "length_zone_median", "length_zone_sd")
    distances <- c(TH, n_zones, dist_zone_min, dist_zone_max, dist_zone_mean,
        dist_zone_median, dist_zone_sd)
    names(distances) <- c("TH", "n_zones", "dist_zone_min", "dist_zone_max",
        "dist_zone_mean", "dist_zone_median", "dist_zone_sd")

    return(list(n_zones = n_zones, n_bases = n_bases, lengths = lengths,
        distances = distances, TH = TH, chr = accumulation$chr,
        w = accumulation$w, acctype = accumulation$acctype))

}
