##Finds transcription factor dense DNA zones for different accumulation
##threshold values.
##input:
##    accumulation: a list of four elements, as the output of the accumulation
##    function.
##    threshold_step: an integer, the step used to calculate the threshold
##    values.
##output: a list of seven elements, and, for each threshold value, a ".bed"
##file with the chromosome and genomic coordinates of the dense zones found:
##    zones_count: a data frame containing the considered threshold values and
##    the number of dense zones obtained with each of the threshold values.
##    bases_count: a data frame containing the considered threshold values and
##    the total number of bases belonging to the dense zones obtained with
##    each of the threshold values.
##    lengths: a data frame containing the considered threshold values and
##    min, max, mean, median and standard deviation of the dense zone lengths
##    obtained with each of the considered threshold values.
##    distances: a dataframe containing the considered threshold values and
##    min, max, mean, median and standard deviation of the distances between
##    adjacent dense zones obtained with each of the threshold values.
##    acctype: a string with the accumulation type used.
##    chr: a string with the chromosome name associated with the accumulation
##    vector used.
##    w: an integer with half-width of the window used to calculate the
##    accumulation vector.




dense_zones <- function(accumulation, threshold_step) {

    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")
    if (class(accumulation$accvector) != "Rle")
        stop("'accvector' element of the list 'accumulation' must be of class
        'Rle'.")
    if (!is.numeric(threshold_step))
        stop("'threshold_step' must be an object of type 'numeric'.")
    if (threshold_step < 1)
        stop("'threshold_step' must be > 0.")
    ## initializing vectors
    n_zones <- vector()
    n_bases <- vector()
    length_zone_min <- vector()
    length_zone_max <- vector()
    length_zone_mean <- vector()
    length_zone_median <- vector()
    length_zone_sd <- vector()
    dist_zone_min <- vector()
    dist_zone_max <- vector()
    dist_zone_mean <- vector()
    dist_zone_median <- vector()
    dist_zone_sd <- vector()
    acc_tf <- accumulation$accvector
    threshold <- seq(1, max(acc_tf), threshold_step)


    for (i in seq_along(threshold)) {

        zones <- slice(acc_tf, i)

        ## counting the number of dense zones that are formed by contiguos bases
        ## given a threshold, and the number of bases belonging to the zones.
        n_bases[i] <- sum(width(zones))
        n_zones[i] <- length(zones)
        ## writing on files chromosome and positions of starting and ending
        ## points
        write.table(data.frame(rep(accumulation$chr, length(zones)),
            start(zones) - 1, end(zones)), file = paste(accumulation$acctype,
            "_acc_w_", accumulation$w, "_", accumulation$chr,
            "_dense_zones_th_", threshold[i], ".bed", sep = ""),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

        ## finding elements of length dataframe

        length_zone_min[i] <- min(width(zones))
        length_zone_max[i] <- max(width(zones))
        length_zone_mean[i] <- mean(width(zones))
        length_zone_median[i] <- median(width(zones))
        length_zone_sd[i] <- sd(width(zones))

        ## calculating distances of dense zones
        if (length(zones) > 1) {

            dist_zone <- start(zones)[-1] - end(zones[-length(end(zones))])
            dist_zone_min[i] <- min(dist_zone)
            dist_zone_max[i] <- max(dist_zone)
            dist_zone_mean[i] <- mean(dist_zone)
            dist_zone_median[i] <- median(dist_zone)
            dist_zone_sd[i] <- sd(dist_zone)
        } else {
            dist_zone <- NA
            dist_zone_min[i] <- NA
            dist_zone_max[i] <- NA
            dist_zone_mean[i] <- NA
            dist_zone_median[i] <- NA
            dist_zone_sd[i] <- NA


        }
    }

    ## result dataframes creation

    zones_count <- data.frame(threshold, n_zones)
    bases_count <- data.frame(threshold, n_bases)
    lengths <- data.frame(threshold, n_zones, length_zone_min, length_zone_max,
        length_zone_mean, length_zone_median, length_zone_sd)
    distances <- data.frame(threshold, n_zones, dist_zone_min, dist_zone_max,
        dist_zone_mean, dist_zone_median, dist_zone_sd)
    return(list(zones_count = zones_count, bases_count = bases_count,
        lengths = lengths, distances = distances, chr = accumulation$chr,
        w = accumulation$w, acctype = accumulation$acctype))

}
