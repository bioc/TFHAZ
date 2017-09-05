\name{high_accumulation_zones}
\alias{high_accumulation_zones}

\title{Finds transcription factor high accumulation DNA zones.
}
\description{
This function finds transcription factor high accumulation DNA zones (TFHAZ). 
Starting from the accumulation vector calculated with the \emph{accumulation}
function, each TFHAZ is formed by contiguous bases with accumulation higher than
the threshold (TH). The threshold is found, considering all the accumulation
values of the vector higher than zero, with the following formula:
\emph{TH = mean(accvector) + 2*std(accvector)}. The function finds also the
number of high accumulation zones, the number of total bases belonging to these
zones, the minimum, maximum, mean, median and standard deviation of these zone
lengths and of the distances between adjacent high accumulation zones.
Furhermore, the function plots, for each chromosome base (x axis), the value of
accumulation (y axis) calculated with the \emph{accumulation} function.
On this graph there are also shown the threshold (with a red line) and, on the x
axis, the bases belonging to the high accumulation zones (with red boxes). 
The plot is saved in a ".png" file.
}
\usage{
high_accumulation_zones(accumulation)
}

\arguments{
\item{accumulation}{
a list of four elements containing: a Rle pbject with accumulation values
(e.g., obtained with the \emph{accumulation} function), the accumulation type,
a chromosome name, and the half-width of the window used for the accumulation
count.
}
}

\value{
A list of eight elements:
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


A ".bed" file with the chromosome and genomic coordinates of the high
accumulation zones found is created.
}

\examples{
# loading data
data("data_man")
# TF_acc_21_w_0 is in the data_man collection of datasets
# to find high accumulation zones
TFHAZ_21_w_0 <- high_accumulation_zones(TF_acc_21_w_0)

}

