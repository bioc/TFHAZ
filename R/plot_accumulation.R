## The function plots the accumulation vector obtained with the accumulation
## function. The function works only if the "accumulation" in input is
##calculated for a single chromosome.
## Input:
##    accumulation: a list of four elements, as the output of the accumulation
##    function.
##    chr: optional argument, needed if the accumulation input is found with chr
##    = "all"; a string or a vector containing strings with the name of the
##    chromosome(s) (e.g., "chr1" or c("chr1", "chr4"))
## Output: a .png file (or more) with the plot(s).

plot_accumulation <- function(accumulation, chr = NULL) {
    if (!is.list(accumulation))
        stop("'accumulation' must be an object of type 'list'.")

    if (class(accumulation$accvector) != "Rle" & class(accumulation$accvector)
        != "SimpleRleList")
        stop("'accvector' element of the list 'accumulation' must be of class
                'Rle' or 'SimpleRleList'.")

    if (accumulation$chr != "all") {

    png(filename = paste(accumulation$acctype, "_w_", accumulation$w, "_",
                        accumulation$chr, ".png", sep = ""))
    plot(accumulation$accvector, type = "l", xlab = "base", ylab =
            paste("# of ", accumulation$acctype, "s", sep = ""))
    dev.off()
    } else {

    if (accumulation$chr == "all" & chr == NULL)
        stop("argument 'chr' is missing")

    for (i in seq_along(chr)) {

        acc <- eval(parse(text = paste("a$accvector$", chr[i], sep = "")))
        png(filename = paste(accumulation$acctype, "_w_", accumulation$w, "_",
        chr[i], ".png", sep = ""))
        plot(acc, type = "l", xlab = "base", ylab =
        paste("# of ", accumulation$acctype,
        "s", sep = ""))
    dev.off()
        }
    }
}
