\name{Ishikawa}
\alias{Ishikawa}
\docType{data}
\title{
Contains genomic regions of transcription factors at the ranges side
and the name of the transcription factors at the metadata side.
}
\description{
\emph{Ishikawa} is a Large GRanges object with 283,009 ranges and 1 metadata column. Each range represents the coordinates of a TF binding region while metadata column indicates the name of the TF.
}
\usage{data("Ishikawa")}
\value{
None, the function is invoked for its side effect.
}
\format{
An object of class "GRanges"
}
\details{
The dataset is obtained from computation of ENCODE ChIP-Seq data of the localization of transcription factor binding sites of the Ishikawa cell line. The data have been processed and extracted with GMQL (GenoMetric Query Language http://www.bioinformatics.deib.polimi.it/ GMQL/).
}
\examples{
# loading dataset
data("Ishikawa")
head(Ishikawa)
}
\keyword{datasets}
