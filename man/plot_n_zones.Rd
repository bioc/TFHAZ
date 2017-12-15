\name{plot_n_zones}
\alias{plot_n_zones}

\title{
Plots the number of dense zones for each threshold value.
}
\description{
For each accumulation threshold value considered, this function shows the number
of dense zones found (e.g., with the \emph{dense_zones} function). It also plots
in red color the point of the graph with maximum slope change (maximum second
derivative).
}
\usage{
plot_n_zones(zones, chr = NULL)
}

\arguments{
\item{zones}{
a list of eight elements, as the output of the \emph{dense_zones} function.
}
\item{chr}{
optional argument, needed if the input was found with chr = "all";
a string or a vector containing strings with the name of the chromosome(s)
(e.g., "chr1" or c("chr1", "chr4")) to be considered
}
}

\value{
None, the function is invoked for its side effect.
}

\examples{
# loading data
data("data_man")
# TF_dense_w_0 is in the data_man collection of datasets
# plot number of zones
plot_n_zones(TF_dense_w_0)
}


