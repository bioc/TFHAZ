## The function is used to plot the number of dense zones and the number of
## bases belonging to these dense zones obtained (all with accumulation
## threshold=1) from different values of w.
## Input:
##    input_list: a list containing multiple outputs of the dense_zones
##    function, obtained for the same accumulation type, same chromosome and
##    different values of w.
##    chr: optional argument, needed if the input is found with chr = "all";
##    a string or a vector containing strings with the name of the chromosome(s)
##    (e.g., "chr1" or c("chr1", "chr4"))

w_analysis <- function(input_list, chr = NULL) {
    if (!is.list(input_list))
        stop("'input_list' must be an object of type 'list'.")

    chromosomes <- vector()
    acc_types <- vector()
    windows <- vector()
    for (i in 1:length(input_list)) {
        chromosomes[i] <- input_list[[i]]$chr
        windows[i] <- input_list[[i]]$w
        acc_types[i] <- input_list[[i]]$acctype
    }
    if (length(unique(windows)) != length(input_list))
        stop("The sizes of the windows are the same. The comparison it's not
            possible")


    ## Check that the input data acc_types are the same
    if (length(unique(acc_types)) != 1)
        stop("The accumulation types related to the input
            arguments are not the same. The comparison it's not possible")

    if (unique(chromosomes) != "all" & length(unique(chromosomes)) == 1) {
        chr <- unique(chromosomes)
    }


    for(k in seq_along(chr)) {

        n_bases <- vector()
        n_zones <- vector()

        for (i in 1:length(input_list)) {

        n_bases[i] <- eval(parse(text = paste("input_list[[i]]$bases_count$",
                                        chr[k],"$n_bases[1]", sep = "")))
        n_zones[i] <- eval(parse(text = paste("input_list[[i]]$zones_count$",
                                        chr[k],"$n_zones[1]", sep = "")))
    }

    ## Check that the input data windows are different




    ## plot(x axis logarithmic)
    plot(windows[order(windows)], n_bases[order(windows)], log = "x",
        type = "l", lty = 2, ylab = "# bases", col = "blue", xaxt = "n", yaxt =
        "n", xlab = "neighborhood half-width (log)", main = paste(chr[k],":
        Number of detected events (DNA zones or bases) \nvs. neighborhood window
        half-width",sep=""))
    axis(1, at = windows)
    axis(2, at = n_bases, col = "blue")
    par(new = TRUE)
    plot(windows[order(windows)], n_zones[order(windows)], log = "x",
        type = "l", ylab = "", xlab = "", col = "red", axes = FALSE)
    axis(4, at = n_zones, col = "red")
    mtext("# zones", side = 4)
    legend("left", legend = c("# bases", "# zones"), col = c("blue", "red"),
        lty = c(2,1), cex = 0.6)

    }
}
