## This function finds transcription factor high accumulation DNA zones (TFHAZ).
## Two different methods for the search of TF high accumulation DNA zones
## are available, with two possible methods to find the threshold value.
##
## input:
##    accumulation: list of four elements containing: a sparse vector with
##    accumulation values (e.g.,obtained with the accumulation function),
##    the accumulation type, a chromosome name, and the half-width of the
##    window used for the accumulation count.
##    method:  a string with the name of the method used to find high
##    accumulation zones: "binding_regions" or "overlaps".
##    data a GRanges object containing coordinates of TF binding regions and
##    their TF name. It is needed in the case of "binding regions" method
##    threshold: a string with the name of the method used to find the threshold
##    value: "std" or "top_perc".
##    k: an integer with the percentage (with the top_perc method) or the
##    number of std deviations (with the std method) to be used in order to 
##    find the threshold.
## output: A list of nine elements:
##    zones: a GRanges object containing the coordinates of the high
##    accumulation zones.
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
## Furthermore, a ".bed" file with the chromosome, genomic coordinates and 
## accumulation values of the found accumulation zones and a ".png" file with 
## the plot of the TFHAZ found along the chromosome (only if the "accumulation" 
## in input is calculated for a single chromosome) can be created.



high_accumulation_zones <- function(accumulation, method = c("overlaps",
                                    "binding_regions"), data, threshold =
                                    c("std", "top_perc"), k, writeBed =
                                        FALSE, plotZones = FALSE)
        {

    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")
    if (class(accumulation$accvector) != "Rle" & class(accumulation$accvector)
        != "SimpleRleList")
        stop("'accvector' element of the list 'accumulation' must be of class
            'Rle'.")
    if (!is.character(method))
        stop("'method' must be an object of type 'character'.")
    if (!is.character(threshold))
        stop("'threshold' must be an object of type 'character'.")
    ## input accumulation vector
    acc_tf <- accumulation$accvector
    if(accumulation$chr != "all") {
        acc_regions <- GRanges(seqnames = rep(accumulation$chr,
            length(ranges(acc_tf))), ranges = ranges(acc_tf), score =
            runValue(acc_tf))
    } else {
    acc_regions <- as(acc_tf,"GRanges") }


    if (method == "binding_regions") {
        if (missing(data))
            stop("argument data is missing")
        if (accumulation$w != 0)
            stop("in binding regions method w must be equal to 0")
        acc_reg <- reduce(acc_regions[score(acc_regions) > 0])
        tf <- unique(elementMetadata(data))[,1]
        count_list <- list()
        for(i in seq_along(tf)) {
            data_tf <- data[elementMetadata(data)[,1] == tf[i]]
            count_list[[i]] <- countOverlaps(acc_reg, data_tf)
            count_list[[i]][count_list[[i]] > 1] <- 1
        }
        count_list <- Reduce("+", count_list)
        acc_regions <- GRanges(acc_reg, score = count_list)
    }

    if (method == "overlaps") {
        acc_regions <- acc_regions[score(acc_regions) > 0]
    }

    ## finding the threshold
    if (threshold == "top_perc") {
        if(missing(perc))
            stop("argument perc is missing")
        if (!is.numeric(perc))
            stop("'perc' must be an object of type 'numeric'.")
        if (perc < 0 || perc >= 100)
            stop("'perc' must be >= 0 and <= 100.")
        n_regions <- length(acc_regions)
        ul <- unique(score(acc_regions))
        max_index <- max(ul)
        ul <- sort(ul,decreasing = FALSE)
        up_perc <- n_regions * (1 - (perc * 0.01))
        sum_ACC_count=vector()
        ACC_count=vector()
        for(i in seq_along(ul)){
            data_ACC <- acc_regions[score(acc_regions) == ul[i]]
            ACC_count[i] <- length(data_ACC)
            if (i == 1) {
                sum_ACC_count[i] <- ACC_count[i]
            } else {
                sum_ACC_count[i]=sum_ACC_count[i-1]+ACC_count[i]
            }
        }
        TH <- ul[which(sum_ACC_count > up_perc)[1]]
        
    }

    if (threshold == "std") {
      if(missing(perc))
        perc <- 2
      acc <- score(acc_regions) * width(acc_regions)
      mean_acc <- sum(acc) / sum(width(acc_regions))
      std_reg = width(acc_regions) * (score(acc_regions) - mean_acc)^2
      std_acc <- sqrt(sum(std_reg) / (sum(width(acc_regions)) - 1))
      TH <- ceiling(mean_acc + perc * std_acc)
      ul <- unique(score(acc_regions))
      max_index <- max(ul)
    }

    ## finding high accumulation zones
    zones_old <- reduce(acc_regions[score(acc_regions) >= TH])
    acc_regions <- acc_regions
    zones <- reduceKeepAttr(acc_regions[score(acc_regions) >= TH], 
                             keep.names = TRUE, with.revmap = TRUE)
    indexes <- as.numeric(zones@elementMetadata@listData[["revmap"]]@partitioning@end)
    acc_val <- score(acc_regions)[indexes]    
    ## finding the number of bases belonging to the zones
    n_bases <- sum(width(zones))
    ## finding the number of zones
    n_zones <- length(zones)




    if (accumulation$chr != "all" & plotZones == TRUE) {
        ## finding the bases belonging to the zones
        bases <- Rle()
        for (i in seq_along(zones)) {
            bases <- append(bases, start(zones)[i]:end(zones)[i])
        }

        ## plot
        png(filename = paste("high_accumulation_zones_TH_", round(TH, digits =
                            1), "_", accumulation$acctype, "_acc_w_",
                            accumulation$w, "_", accumulation$chr, ".png",
                            sep = ""))
        plot(accumulation$accvector, type = "l", xlab = "base", ylab = paste
            ("# of ", accumulation$acctype, "s", sep = ""))
        points(bases, rep(0, length(bases)), col = "red", pch = 15)
        abline(h = TH, col = "red")
        legend("topleft", legend = c("threshold", "HOT\nzones"),
                col = "red", pch = c(NA, 15), lty = c(1, NA), bty='n')
        dev.off()
    }

    if (writeBed == TRUE) {
    ## writing on files chromosome and positions of starting and ending points
    write.table(data.frame(seqnames(zones), start(zones) - 1, end(zones), acc_val),
                file = paste(accumulation$acctype, "_acc_w_", accumulation$w,
                "_", accumulation$chr, "_dense_zones_th_", round(TH, digits = 1)
                , ".bed", sep = ""), row.names = FALSE, col.names = FALSE,
                quote = FALSE, sep = "\t") }

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

    return(list(zones = zones, n_zones = n_zones, n_bases = n_bases, lengths =
                lengths, distances = distances, TH = TH, chr = accumulation$chr,
                w = accumulation$w, max_accumulation_index = max_index, acctype = accumulation$acctype))

}
