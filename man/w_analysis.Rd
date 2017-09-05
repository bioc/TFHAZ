\name{w_analysis}
\alias{w_analysis}

\title{
Shows the number of dense zones obtained with different values of base
neighborhood window width.
}
\description{
This function is used to plot the number of dense zones and the total number of
bases belonging to these dense zones present in a set of inputs, obtained (all
with accumulation threshold=1) using the \emph{dense_zones} function, for the
same accumulation type, same chromosome, and different values of \emph{w}
half-width of the window defining the neighborhood of each base. The plot (with
x axis logarithmic-scale) is generated only if all input data refer to the same
chromosome and accumulation type, otherwise an error message appears. Beside
helping in comparing results obtained with different values of \emph{w}, this
function supports finding the best value of \emph{w} to be used in further
analysis.
}
\usage{
w_analysis(input_list)
}

\arguments{
\item{input_list}{
a list containing multiple outputs of the \emph{dense_zones} function, obtained
for the same accumulation type, same chromosome and different values of
\emph{w}.
}
}
\value{
None, the function is invoked for its side effect.
}
\examples{
# loading data
data("data_man")
# l is a list with dense zone (e.g., TF_dense_21_w_10, TF_dense_21_w_100,
# TF_dense_21_w_1000 and TF_dense_21_w_10000) objects present in data_man
# collection of datasets
l <- list(TF_dense_21_w_10, TF_dense_21_w_100, TF_dense_21_w_1000,
TF_dense_21_w_10000)
# plot
w_analysis(l)
}


