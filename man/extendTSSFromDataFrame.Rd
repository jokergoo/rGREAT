\name{extendTSSFromDataFrame}
\alias{extendTSSFromDataFrame}
\title{
Extend TSS
}
\description{
Extend TSS
}
\usage{
extendTSSFromDataFrame(df, seqlengths, genome = NULL,
    strand = NULL, gene_id = NULL,
    gene_id_type = NULL, verbose = great_opt$verbose, ...)
}
\arguments{

  \item{df}{A bed-like data frame where the first three columns should be chromosomes, start positions, end positions. It does not matter whether regions correspond to genes or TSS. }
  \item{seqlengths}{A named vector of chromosome lengths.}
  \item{genome}{UCSC genome can be set here, then \code{seqlengths} will be automatically retrieved from UCSC server.}
  \item{strand}{The strand information can be provided in \code{df} as a column named "strand" or as a column with "+"/"-"/"*", or the strand information can be provided as a vector and be assigined to this argument.}
  \item{gene_id}{The gene ID information can be provided in \code{df} as a column named "gene_id", or it can be provided as a vector and be assigned to this argument.}
  \item{gene_id_type}{Gene ID types in \code{df}. You need to set this argument if you use built-in gene sets in \code{\link{great}} so that genes can be correctly mapped. The value can only be one of "SYMBOL", "ENTREZ", "ENSEMBL" and "REFSEQ".}
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
