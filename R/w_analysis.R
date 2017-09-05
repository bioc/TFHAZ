## The function is used to plot the number of dense zones and the number of
## bases belonging to these dense zones obtained (all with accumulation
## threshold=1) from different values of w.
## Input: a list containing multiple outputs of the dense_zones function,
## obtained for the same accumulation type, same chromosome and different values
## of w.
## Output: The two plot, the second one with logarithmic-scale x axis.
w_analysis <- function(input_list) {
    if (!is.list(input_list))
        stop("'input_list' must be an object of type 'list'.")
    ## Initialize vectors
    windows <- vector()
    chromosomes <- vector()
    acc_types <- vector()
    n_bases <- vector()
    n_zones <- vector()

    ## vectors with input data
    for (i in 1:length(input_list)) {
        windows[i] <- input_list[[i]]$w
        chromosomes[i] <- input_list[[i]]$chr
        acc_types[i] <- input_list[[i]]$acctype
        n_bases[i] <- input_list[[i]]$bases$n_bases[1]
        n_zones[i] <- input_list[[i]]$zones$n_zones[1]
    }

    ## Check that the input data windows are different
    if (length(unique(windows)) != length(input_list)) 
        stop("The sizes of the windows are the same. The comparison it's not
            possible")
    

    ## Check that the input data chromosomes or data acc_types are the same
    if (length(unique(chromosomes)) != 1 || length(unique(acc_types)) != 1)
        stop("The chromosomes or the accumulation types related to the input
            arguments are not the same. The comparison it's not possible")
    


    ## plot(x axis logarithmic)
    plot(windows[order(windows)], n_bases[order(windows)], log = "x",
        type = "l", lty = 2, ylab = "# bases", col = "blue", xaxt = "n", yaxt =
        "n", xlab = "neighborhood half-width (log)", main = "Number of detected
        events (DNA zones or bases) \nvs. neighborhood window half-width")
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
