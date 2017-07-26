\name{plot_n_zones}
\alias{plot_n_zones}

\title{
Plots the number of dense zones for each threshold value.
}
\description{
For each accumulation threshold value considered, this function shows the number of dense zones found (e.g., with the \emph{dense_zones} function). It also plots in red color the point of the graph with maximum slope change (maximum second derivative).
}
\usage{
plot_n_zones(zones)
}

\arguments{
  \item{zones}{
a list of seven elements, as the output of the \emph{dense_zones} function.
}
}

\value{
None, the function is invoked for its side effect.
}

\examples{
# loading data
data("data_man")
# TF_dense_21_w_0 is in the data_man collection of datasets
# plot number of zones
plot_n_zones(TF_dense_21_w_0)
}


