\name{n_zones_PCA}
\alias{n_zones_PCA}

\title{
Principal Component Analysis of the number of dense zones obtained with the three methods of accumulation (TF, region, base).
}
\description{
This function performs the Principal Component Analysis (PCA) of the number of dense zones obtained by varying the threshold on accumulation values obtained with the three methods of accumulation (\emph{TF}, \emph{region}, \emph{base}). Before performing the PCA, the number of dense zone values are scaled with the \emph{scale} R function. This function works only if the number of different threshold values used to find the dense zones with the \emph{dense_zones} function is the same for all the three accumulation types, while the threshold values can be different.
}
\usage{
n_zones_PCA(TF_zones, region_zones, base_zones)
}

\arguments{
  \item{TF_zones}{
a list with the results of the \emph{dense_zones} function using the TF accumulation method and varying the thresholds on the considered accumulation values.
}
  \item{region_zones}{
a list with the results of the \emph{dense_zones} function using the region accumulation method and varying the thresholds on the considered accumulation values.
}
  \item{base_zones}{
a list with the results of the \emph{dense_zones} function using the base accumulation method and varying the thresholds on the considered accumulation values.
}
}
\value{
A list with a summary containing the standard deviation on each principal component, the proportion of variance explained by each principal component, the cumulative proportion of variance described by each principal component, and the loadings of each principal component. In addition, a plot with the variances of the principal components; a plot with the cumulate variances of the principal components and a plot with the loadings of the three principal components.


}

\examples{
# loading data
data("data_man")
# TF_dense_21_w_10, reg_dense_21_w_10, and base_dense_21_w_10 are in the data_man
# collection of datasets
# PCA
n_zones_PCA(TF_dense_21_w_10, reg_dense_21_w_10, base_dense_21_w_10)
}


