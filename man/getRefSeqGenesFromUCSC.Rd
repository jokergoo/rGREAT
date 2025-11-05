\name{getRefSeqGenesFromUCSC}
\alias{getRefSeqGenesFromUCSC}
\title{
Get RefSeq genes from UCSC
}
\description{
Get RefSeq genes from UCSC
}
\usage{
getRefSeqGenesFromUCSC(genome, subset = c("RefSeqSelect", "RefSeqCurated"))
}
\arguments{

  \item{genome}{UCSC genome, such as "hg19".}
  \item{subset}{Subset of RefSeq genes. See \url{https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=refSeqComposite} .}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object.
}
\examples{
# There is no example
NULL

}
