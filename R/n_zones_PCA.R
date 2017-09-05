## The function performs the Principal Component Analysis (PCA) of the number of
## dense zones obtained by varying the threshold on accumulation values obtained
## with the three methods of accumulation (TF, region, base).
## Input: 3 lists with the results of the dense_zones function using the TF,
## region and base accumulation method (in this order) and varying the
## thresholds on the considered accumulation values.
## Output: A list with a summary containing the standard deviation on each
## principal component, the proportion of variance explained by each principal
## component, the cumulative proportion of variance described by each principal
## component, and the loadings of each principal component. In addition, a plot
## with the variances of the principal components; a plot with the cumulate
## variances of the principal components; three plots with the loadings of the
## three principal components; a plot with the vector obtained, for each
## threshold value, as linear combination of the number of dense zones
## (after scaling) for each accumulation type, using as linear combination
## coefficients the loadings of the first principal component for the three
## accumulation types.

n_zones_PCA <- function(TF_zones, region_zones, base_zones) {
    if (!is.list(TF_zones))
        stop("'TF_zones' must be an object of type 'list'.")
    if (!is.list(region_zones))
        stop("'region_zones' must be an object of type 'list'.")
    if (!is.list(base_zones))
        stop("'base_zones' must be an object of type 'list'.")
    if (TF_zones$acctype != "TF")
        stop("TF_zones must contain dense zones from TF accumulation type.")
    if (region_zones$acctype != "region")
    stop("region_zones must contain dense zones from region accumulation type.")
    if (base_zones$acctype != "base")
    stop("base_zones must contain dense zones from base accumulation type.")
    if (length(TF_zones$zones_count$n_zones) !=
        length(region_zones$zones_count$n_zones) ||
        length(TF_zones$zones_count$n_zones) !=
        length(base_zones$zones_count$n_zones) ||
        length(region_zones$zones$n_zones) != length(base_zones$zones$n_zones))
        stop("The lengths of the three input dataframe are different,
            comparison not possible")

    ## dataframe with the number of zones for the three inputs
    n_zones <- data.frame(TF_zones$zones_count$n_zones,
        region_zones$zones_count$n_zones, base_zones$zones_count$n_zones)
    ## scaled values of n_zones
    zone <- data.frame(scale(n_zones))

    ## PCA
    pca <- prcomp(zone)
    summary(pca)
    pca$rotation[, 1:2]  # loadings
    ## barplot with variances
    barplot(pca$sdev^2, main = "Variances of the principal components",
        ylab = "Variances", names.arg = c("Comp1", "Comp2", "Comp3"))
    ## plot with cumulate variances
    plot(cumsum(pca$sdev^2/sum(pca$sdev^2)), type = "l", main =
        "Cumulate variances of the principal components", ylab = "Contribute to
        total variance", xlab = "Number of components", xaxt = "n")
    axis(1, c(1, 2, 3))

    ## plots with loadings
    par(mfrow = c(3, 1))
    barplot(pca$rotation[, 1], ylim = c(-1, 1), names.arg = c("TF", "region",
        "base"), ylab = "Comp1")
    barplot(pca$rotation[, 2], ylim = c(-1, 1), names.arg = c("TF", "region",
        "base"), ylab = "Comp2")
    barplot(pca$rotation[, 3], ylim = c(-1, 1), names.arg = c("TF", "region",
        "base"), ylab = "Comp3")
    mtext("Loadings of the principal components", outer = TRUE)
    return(list(summary = summary(pca), loadings = pca$rotation))
}


