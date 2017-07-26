\name{accumulation}
\alias{accumulation}
\title{
Creates a vector with accumulation counts of transcription factors for each chromosome base.
}
\description{
From a dataset with transcription factor (TF) binding regions, this function creates a vector (\emph{accvector}) in which, for each chromosome base, the accumulation of the TFs present in the input dataset is calculated. Three types of accumulation are considered: \emph{TF accumulation}, \emph{region accumulation} and \emph{base accumulation}.
TF accumulation: for each base, it is the number of different TFs present in the neighborhood of the considered base. The neighborhood is defined by a window with half-width \emph{w} centered on the considered base.
Region accumulation: for each base, it is the number of regions containing TFs in the neighborhood of the considered base. If in the neighborhood of a base there are two input binding regions of the same TF, the accumulation value in that base is equal to 2 (differently from the TF accumulation, whose value in the same case is equal to 1).
Base accumulation: for each base, it is the total number of bases belonging to input regions containing TFs in the neighborhood of the considered base.
With w=0, a single base approach is applied (no base neighborhood is considered). In this case, if in the input dataset overlapping regions for the same TF and chromosome do not exist, the results of TF, region and base accumulation are equal.



}
\usage{
accumulation(data, acctype, chr, w)
}

\arguments{
  \item{data}{
a GRanges object containing coordinates of TF binding regions and their TF name.
}
  \item{acctype}{
a string with the name of the accumulation type: "TF", "region", "base".
}
  \item{chr}{
a string with the name of the chromosome (e.g., "chr1").
}
  \item{w}{
an integer, half-width of the window that defines the neighborhood of each base.
}
}

\value{
A list of four elements:
\item{accvector}{a sparse vector containing the accumulation for each base of the selected chromosome.}
\item{acctype}{a string with the accumulation type used.}
\item{chr}{a string with the chromosome name associated with the accumulation vector.}
\item{w}{an integer with the half-width of the window used to calculate the accumulation vector.}
}
\examples{
\donttest{
# loading dataset
data("Ishikawa")
# to calculate TF accumulation for the chromosome 21 for w=0
TF_acc_21_w_0 <- accumulation(Ishikawa,"TF","chr21",0)}

}


