\name{high_accumulation_zones}
\alias{high_accumulation_zones}

\title{Finds transcription factor high accumulation DNA zones.
}
\description{
This function finds transcription factor high accumulation DNA zones (TFHAZ).
Starting from the accumulation vector calculated with the \emph{accumulation}
function, two different methods for the search of TF high accumulation DNA zones
are available. The \emph{binding regions} method is based on the identification
of DNA regions with presence of TF binding (at least one TF) from which those
with a high number of different TFs (above the threshold) are selected. This
method works only if the accumulation vector is found with \emph{w}=0. The
\emph{overlaps} method is the method used also in \emph{dense_zones} function.
It uses a single base local approach, identifying DNA bases, that form the dense
zones, in which there is high overlap of TFs. For the \emph{binding regions}
method the high accumulation zones are the accumulation regions with values
higher or equal to the threshold, while in the \emph{overlaps} these zones are
defined as sets of contiguous bases with accumulation value higher or equal to
the considered threshold. The threshold value is found considering two methods.
The \emph{std} method considers all and only the bases of the accumulation
vector (\emph{accvector}) with values higher than zero, and the threshold is
found with the following formula:
\emph{TH = mean(accvector) + 2*std(accvector)}. The \emph{top_perc} method
considers the accumulation regions and selects those in the top x percentage,
with x chosen by the user through the \emph{perc} argument. The function finds
also the number of high
accumulation zones, the number of total bases belonging to these zones, the
minimum, maximum, mean, median and standard deviation of these zone
lengths and of the distances between adjacent high accumulation zones. In the
case of \emph{binding regions} method, it is needed to include the \emph{data}
input argument, that is the GRanges object used in the \emph{accumulation}
function.
Furhermore, in the case of single chromosome accumulation vector, the function
can plots, for
each chromosome base (x axis), the value of accumulation (y axis) calculated
with the \emph{accumulation} function. On this graph there are also shown the
threshold (with a red line) and, on the x axis, the bases belonging to the high
accumulation zones (with red boxes). The plot can be saved in a ".png" file.
}
\usage{
high_accumulation_zones(accumulation, method = c("overlaps", "binding_regions"),
data, threshold = c("std","top_perc"), perc, writeBed = FALSE, plotZones =
FALSE)
}

\arguments{
\item{accumulation}{
a list of four elements containing: a Rle object (or SimpleRleList) with
accumulation values
(e.g., obtained with the \emph{accumulation} function), the accumulation type,
a chromosome name, and the half-width of the window used for the accumulation
count.
}
\item{method}{
a string with the name of the method used to find high accumulation zones:
"binding_regions" or "overlaps".
}
\item{data}{
a GRanges object containing coordinates of TF binding regions and their TF name.
It is needed in the case of \emph{binding regions} method.
}
\item{threshold}{
a string with the name of the method used to find the threshold value:
"std" or "top_perc".
}
\item{perc}{
an integer with the percentage value to be used in order to find the threshold
with the \emph{top_perc} method.
}
\item{writeBed}{
When set to TRUE, for each threshold value a ".bed"
file with the chromosome and genomic coordinates of the dense zones found is
created.
}
\item{plotZones}{
When set to TRUE, and the "accumulation" in input is calculated for a single
chromosome, a ".png" file with the plot of the high accumulation zones on the
accumulation vector is created.
}
}

\value{
A list of nine elements:
\item{zones}{a GRanges object containing the coordinates of the high
accumulation zones.}
\item{n_zones}{an integer containing the number of high accumulation zones
obtained.}
\item{n_bases}{an integer containing the total number of bases belonging to the
high accumulation zones obtained.}
\item{lengths}{a vector containing the considered threshold value and min, max,
mean, median and standard deviation of the high accumulation zone lengths
obtained.}
\item{distances}{a vector containing the considered threshold value and  min,
max, mean, median and standard deviation of the distances between adjacent high
accumulation zones obtained.}
\item{TH}{a number with the threshold value found.}
\item{acctype}{a string with the accumulation type used.}
\item{chr}{a string with the chromosome name associated with the accumulation
vector used.}
\item{w}{an integer with half-width of the window used to calculate the
accumulation vector.}

}

\examples{
# loading dataset
data("Ishikawa")
# TF_acc_w_0 is in the data_man collection of datasets
# to find high accumulation zones
TFHAZ_w_0 <- high_accumulation_zones(TF_acc_w_0, method = "overlaps",
threshold = "std")

}

