\name{plot_accumulation}
\alias{plot_accumulation}

\title{
Plots the accumulation vector obtained with the accumulation function.
}
\description{
For each chromosome base (x axis), this function plots the value of accumulation
(y axis) calculated with the \emph{accumulation} function. The plot is saved in
a ".png" file.
}
\usage{
plot_accumulation(accumulation)
}

\arguments{
\item{accumulation}{
a list of four elements, as the output of the \emph{accumulation} function,
containing: a sparse vector with accumulation values, the accumulation type, a
chromosome name, and the half-width of the window used for the accumulation
count.
}
}

\value{
None, the function is invoked for its side effect.
}

\examples{
\donttest{
# loading dataset
data("data_man")
# TF_acc_21_w_0 is in the data_man collection of datasets
# plot accumulation vector
plot_accumulation(TF_acc_21_w_0)}
}

