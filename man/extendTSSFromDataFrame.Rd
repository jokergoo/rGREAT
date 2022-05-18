\name{extendTSSFromDataFrame}
\alias{extendTSSFromDataFrame}
\title{
Extend TSS
}
\description{
Extend TSS
}
\usage{
extendTSSFromDataFrame(df, seqlengths, strand = NULL, gene_id = NULL,
    verbose = great_opt$verbose, ...)
}
\arguments{

  \item{df}{A bed-like data frame where the first three columns should be chromosomes, start positions, end positions. It does not matter whether regions correspond to genes or TSS. }
  \item{seqlengths}{A named vector of chromosome lengths.}
  \item{strand}{The strand information can be provided in \code{df} as a column named "strand" or as a column with "+"/"-"/"*", or the strand information can be provided as a vector and be assigined to this argument.}
  \item{gene_id}{The gene ID information can be provided in \code{df} as a column named "gene_id", or it can be provided as a vector and be assigned to this argument.}
  \item{verbose}{Whether to print messages.}
  \item{...}{All pass to \code{\link{extendTSS}}.}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with one meta column 'gene_id'.
}
\examples{
# There is no example
NULL

}
