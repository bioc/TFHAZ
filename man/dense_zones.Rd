\name{dense_zones}
\alias{dense_zones}

\title{
Finds transcription factor dense DNA zones for different accumulation threshold
values.
}
\description{
For each accumulation threshold value defined, this function finds transcription
factor (TF) dense DNA zones (regions). Starting from the accumulation vector
calculated with the \emph{accumulation} function, each dense zone is formed by
contiguous bases with accumulation equal or higher than the threshold.
For each defined threshold value, the function finds also the number of dense
zones, the number of total bases belonging to the dense zones, the minimum,
maximum, mean, median and standard deviation of the dense zone lengths and of
the distances between adjacent dense zones.
}
\usage{
dense_zones(accumulation, threshold_step, chr = NULL, writeBed = FALSE)
}

\arguments{
\item{accumulation}{
a list of four elements containing: a Rle (or SimpleRleList) object with
accumulation values
(e.g., obtained with the \emph{accumulation} function), the accumulation type,
a chromosome name, and the half-width of the window used for the accumulation
count.
}
\item{threshold_step}{
an integer, the step used to calculate the threshold values. These values vary
from 1 to the maximum accumulation value in the considered accumulation vector
(e.g., found with the \emph{accumulation} function). In the case of using
accumulation values of base accumulation type and w different from zero (e.g.,
w=1000), choosing a step quite different from 1 is suggested because the maximum
accumulation value is usually high.
}
\item{chr}{
optional argument, a string with a cromosome name. It is needed to apply the
function only to a single cromosome present in the accumulation input. If
chr = "all" (default value) the function operates on all the
chromosomes present in the input.
}
\item{writeBed}{
When set to TRUE, for each threshold value (and for each chromosome) a ".bed"
file with the chromosome and genomic coordinates of the dense zones found is
created.
}
}
\value{
A list of eight elements:
\item{zones}{a list with "GRanges" objects with the dense zones found for each
chromosome and threshold value considered.
}
\item{zones_count}{a list with a data frame for each chromosome considered,
containing the considered threshold values and the number of dense zones
obtained with each of the threshold values.}
\item{bases_count}{a list with a data frame for each chromosome considered,
containing the considered threshold values and the total number of bases
belonging to the dense zones obtained with each of the threshold values.}
\item{lengths}{a list with a data frame for each chromosome considered,
containing the considered threshold values and min, max, mean, median and
standard deviation of the dense zone lengths obtained with each of the
considered threshold values.}
\item{distances}{a list with a data frame for each chromosome considered,
containing the considered threshold values and min, max, mean, median and
standard deviation of the distances between adjacent dense zones obtained with
each of the threshold values.}
\item{acctype}{a string with the accumulation type used.}
\item{chr}{a string with the chromosome name associated with the output zones.}
\item{w}{an integer with half-width of the window used to calculate the
accumulation vector.}

When writeBed is set to TRUE, for each threshold value (and for each chromosome)
a ".bed" file with the chromosome and genomic coordinates of the dense zones
found is created.

}

\examples{
# loading data
data("data_man")
# TF_acc_w_0 is in the data_man collection of datasets
# to find dense zones, with threshold step equal to 1
TF_dense_w_0 <- dense_zones(TF_acc_w_0, 1)
}


