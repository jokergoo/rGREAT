\name{getGapFromUCSC}
\alias{getGapFromUCSC}
\title{
Get gap regions from UCSC
}
\description{
Get gap regions from UCSC
}
\usage{
getGapFromUCSC(genome)
}
\arguments{

  \item{genome}{UCSC genome, such as "hg19".}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object.
}
\examples{
getGapFromUCSC("hg19")
}
