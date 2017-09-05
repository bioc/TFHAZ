## The function plots the accumulation vector obtained with the accumulation
## function.
## Input: a list of four elements, as the output of the accumulation function.
## Output: a .png file with the plot.

plot_accumulation <- function(accumulation) {
    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")
    if (class(accumulation$accvector) != "Rle")
        stop("'accvector' element of the list 'accumulation' must be of class
        'Rle'.")
    png(filename = paste(accumulation$acctype, "_w_", accumulation$w, "_",
        accumulation$chr, ".png", sep = ""))
    plot(accumulation$accvector, type = "l", xlab = "base", ylab =
        paste("# of ", accumulation$acctype,
        "s", sep = ""))
    dev.off()
}
