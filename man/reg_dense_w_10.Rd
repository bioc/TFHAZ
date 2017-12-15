\name{reg_dense_w_10}
\alias{reg_dense_w_10}
\docType{data}
\title{
Contains an output of the dense_zones function.
}
\description{
\emph{reg_dense_w_10} is a list of 8 elements (\emph{zones}, \emph{zones_count},
\emph{bases_count}, \emph{lengths}, \emph{distances}, \emph{acctype},
\emph{chr}, \emph{w}). It is the output of the \emph{dense_zones} function
applied to the accumulation vector found with
w=10, chr="all", acctype="reg". \emph{reg_dense_w_10} is included in the
\emph{data_man} collection.
}
\usage{data("data_man")}
\value{
None, the function is invoked for its side effect.
}
\format{
An object of class list.
}
\examples{
# loading data
data("data_man")
# reg_dense_w_10 is in the data_man collection of datasets
head(reg_dense_w_10)
}
\keyword{datasets}
