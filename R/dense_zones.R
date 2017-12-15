##Finds transcription factor dense DNA zones for different accumulation
##threshold values. The function works only if the "accumulation" in input is
##calculated for a single chromosome.
##input:
##    accumulation: a list of four elements, as the output of the accumulation
##    function.
##    threshold_step: an integer, the step used to calculate the threshold
##    values.
##    chr: optional argument, a string with a cromosome name. It is needed to
##    apply the function only to a single cromosome present in the accumulation
##    input found with chr = "all". If chr = "all" (default value) the function
##    operates for all the chromosome present in the input.
##    writeBed: When set to TRUE, for each threshold value (and for each
##    chromosome) a ".bed" file with the chromosome and genomic coordinates of
##    the dense zones found is created.
##output: a list of eight elements:
##    zones: a list with "GRanges" objects with the dense zones found for each
##    chromosome and threshold value considered.
##    zones_count: a list with a data frame, for each chromosome considered,
##    containing the considered threshold values and the number of dense zones
##    obtained with each of the threshold values.
##    bases_count: a list with a data frame, for each chromosome considered,
##    containing the considered threshold values and the total number of bases
##    belonging to the dense zones obtained with each of the threshold values.
##    lengths: a list with a data frame, for each chromosome considered,
##    containing the considered threshold values and min, max, mean, median and
##    standard deviation of the dense zone lengths obtained with each of the
##    considered threshold values.
##    distances: a list with a data frame, for each chromosome considered,
##    containing the considered threshold values and min, max, mean, median and
##    standard deviation of the distances between adjacent dense zones obtained
##    with each of the threshold values.
##    acctype: a string with the accumulation type used.
##    chr: a string with the chromosome name associated with the accumulation
##    vector used.
##    w: an integer with half-width of the window used to calculate the
##    accumulation vector.
##    When writeBed is set to TRUE, for each threshold value (and for each
##    chromosome) a ".bed" file with the chromosome and genomic coordinates of
##    the dense zones found is created.



dense_zones <- function(accumulation, threshold_step, chr = NULL, writeBed =
                            FALSE) {

    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")

    if (class(accumulation$accvector) != "Rle" & class(accumulation$accvector)
        != "SimpleRleList")
        stop("'accvector' element of the list 'accumulation' must be of class
            'Rle' or 'SimpleRleList'.")
    if (!is.numeric(threshold_step))
        stop("'threshold_step' must be an object of type 'numeric'.")
    if (threshold_step < 1)
        stop("'threshold_step' must be > 0.")
    ## initializing vectors
    d_zones <- list()
    zones_count <- list()
    bases_count <- list()
    lengths <- list()
    distances <- list()

    acc_tf <- accumulation$accvector
    if(accumulation$chr != "all") {
        acc_regions <- GRanges(seqnames = rep(accumulation$chr,
                        length(ranges(acc_tf))), ranges = ranges(acc_tf),
                        score = runValue(acc_tf))
        c <- accumulation$chr
    } else {

        acc_regions <- as(acc_tf,"GRanges")

        if(!is.null(chr)) {
            c <- chr
            acc_regions <- acc_regions[seqnames(acc_regions) == chr]
        } else {
            c <- "all"
        }
    }
    acc_regions <- acc_regions[score(acc_regions) > 0]

    u=unique(seqnames(acc_regions))

    for (k in seq_along(u)) {

        acc_regions_chr <- acc_regions[seqnames(acc_regions) == u[k]]
        threshold <- seq(1, max(score(acc_regions_chr)), threshold_step)
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
        z=list()
        for (i in seq_along(threshold)) {

            zones <- reduce(acc_regions_chr[score(acc_regions_chr) >=
                                                threshold[i]])

            z[[i]]=zones
        ## counting the number of dense zones that are formed by contiguos bases
        ## given a threshold, and the number of bases belonging to the zones.

            n_bases[i] <- sum(width(zones))
            n_zones[i] <- length(zones)

            if (writeBed == TRUE) {
            ## writing on files chromosome and positions of starting and ending
            ## points
                write.table(data.frame(rep(u[k], length(zones)),
                            start(zones) - 1, end(zones)),
                            file = paste(accumulation$acctype, "_acc_w_",
                            accumulation$w, "_", u[k], "_dense_zones_th_",
                            threshold[i], ".bed", sep = ""), row.names = FALSE,
                            col.names = FALSE, quote = FALSE, sep = "\t")

                ## finding elements of length dataframe
            }
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
        zones_count[[k]] <- data.frame(threshold, n_zones)
        bases_count[[k]] <- data.frame(threshold, n_bases)

        lengths[[k]] <- data.frame(threshold, n_zones, length_zone_min,
                                length_zone_max, length_zone_mean,
                                length_zone_median, length_zone_sd)

        distances[[k]] <- data.frame(threshold, n_zones, dist_zone_min,
                                dist_zone_max, dist_zone_mean,
                                dist_zone_median, dist_zone_sd)
        names(z) <- paste("th_", threshold, sep="")
        d_zones[[k]] <- z
    }

    ## result dataframes creation

    names(zones_count) <- u
    names(bases_count) <- u
    names(lengths) <- u
    names(distances) <- u
    names(d_zones) <- u
    return(list(zones = d_zones, zones_count = zones_count, bases_count =
                    bases_count, lengths = lengths, distances = distances,
                chr = c, w = accumulation$w, acctype =
                    accumulation$acctype))

}
