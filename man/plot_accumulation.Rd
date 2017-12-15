\name{plot_accumulation}
\alias{plot_accumulation}

\title{
Plots the accumulation vector obtained with the accumulation function.
}
\description{
For each chromosome base (x axis), this function plots the value of accumulation
(y axis) calculated with the \emph{accumulation} function. The plot is saved in
a ".png" file. If the accumulation input was found with chr = "all", the
chromosomes (one or more) to be considered can be chosen.
}
\usage{
plot_accumulation(accumulation, chr = NULL)
}

\arguments{
\item{accumulation}{
a list of four elements, as the output of the \emph{accumulation} function,
containing: a sparse vector with accumulation values, the accumulation type, a
chromosome name, and the half-width of the window used for the accumulation
count.
}
\item{chr}{
optional argument, needed if the accumulation input was found with chr = "all";
a string or a vector containing strings with the name of the chromosome(s)
(e.g., "chr1" or c("chr1", "chr4")) to be considered.
}
}

\value{
None, the function is invoked for its side effect.
}

\examples{
\donttest{
# loading dataset
data("data_man")
# TF_acc_w_0 is in the data_man collection of datasets
# plot accumulation vector
plot_accumulation(TF_acc_w_0)}
}

