## Creates a vector with accumulation counts of transcription factors for each
## chromosome base.
## input:
## data: a GRanges object containing coordinates of TF binding regions
## and their TF name.
## acctype: a string with the name of the accumulation type: "TF",
## "region", "base".
## chr: a string with the name of the chromosome (e.g., "chr1").
## w: an integer, half-width of the window that defines the neighborhood
## of each base.
## output: A list of four elements.
## accvector: a Rle object containing the accumulation for each base of
## the selected chromosome.
## acctype: a string with the accumulation type used.
## chr: a string with the chromosome name associated with the accumulation
## vector.
## w: an integer with the half-width of the window used to calculate the
## accumulation vector.

accumulation <- function(data, acctype = c("TF", "region", "base"), chr, w) {
    if (class(data) != "GRanges")
        stop("'data' must be an object of class 'GRanges'.")
    if (!is.character(acctype))
        stop("'acctype' must be an object of type 'character'.")
    if (!is.character(chr))
        stop("'chr' must be an object of type 'character'.")
    if (!is.numeric(w))
        stop("'w' must be an object of type 'numeric'.")
    if (w < 0)
        stop("'w' must be >= 0.")
    # dataset with only the selected chromosome
    chrom <- data[seqnames(data) == chr]
    acctype <- match.arg(acctype)

    if (acctype == "TF") {
        # TF accumulation

        # find the number of different TF present in the chromosome
        tf <- unique(elementMetadata(data))[, 1]
        acc_tf <- list()
        for (k in seq_along(tf)) {

            # dataset for a single TF and the selected chromosome
            chrom_tf <- data[seqnames(data) == chr & elementMetadata(data)[, 1]
                == tf[k]]
            # ranges expansion
            start(chrom_tf) <- pmax(1, start(chrom_tf) - w)
            end(chrom_tf) <- end(chrom_tf) + w
            # calculating the coverage
            v_tf <- coverage(chrom_tf, width = max(end(chrom)) + w)
            # max accumulation possible value = 1
            v_tf[v_tf >= 1] <- 1

            acc_tf[[k]] <- v_tf

        }
        # sum of all different TFs
        acc <- Reduce(`+`, acc_tf)
        acc <- eval(parse(text = paste("acc$", chr, sep = "")))
        return(list(accvector = acc, chr = chr, w = w, acctype = acctype))
    }

    if (acctype == "region") {
        # region accumulation
        # ranges expansion
        start(chrom) <- pmax(1, start(chrom) - w)
        end(chrom) <- end(chrom) + w
        # calculating the coverage
        v_r <- coverage(chrom, width = max(end(chrom)))
        acc <- eval(parse(text = paste("v_r$", chr, sep = "")))
        return(list(accvector = acc, chr = chr, w = w, acctype = acctype))
    }

    if (acctype == "base") {
        # base accumulation
        # calculating the coverage
        v_b <- coverage(chrom, width = max(end(chrom)))
        v_b <- eval(parse(text = paste("v_b$", chr, sep = "")))
        # calculating the sum of the bases belonging to the neighborhood
        v_b_1 <- append(v_b, rep(0, w))
        v_b_2 <- append(rev(v_b), rep(0, w))
        acc <- runsum(v_b_1, w + 1) + rev(runsum(v_b_2, w + 1)) - v_b
        return(list(accvector = acc, chr = chr, w = w, acctype = acctype))
    }

}

