\name{extendTSSFromTxDb}
\alias{extendTSSFromTxDb}
\title{
Extend TSS
}
\description{
Extend TSS
}
\usage{
extendTSSFromTxDb(txdb, chromosomes = NULL,
    verbose = great_opt$verbose, ...)
}
\arguments{

  \item{txdb}{Name of "TxDb.*" packages from Bioconductor. All supported TxDb packages are in \code{rGREAT:::BIOC_ANNO_PKGS$txdb}. Note short genome version can also be used here such as "hg19" or "hg19.knownGene".}
  \item{chromosomes}{The original TxDb data contains not only the "normal chromosomes" but also a lot of small contigs. This function can automatically identify those "normal chromosomes" (not guaranteed), but users can always control which chromosomes to use by specify a vector of chromosome names to this argument.}
  \item{verbose}{Whether to print messages.}
  \item{...}{All pass to \code{\link{extendTSS}}.}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with one meta column 'gene_id'.
}
\examples{
if(FALSE) {
extendTSSFromTxDb("TxDb.Hsapiens.UCSC.hg19.knownGene")
extendTSSFromTxDb("hg19")
}
}
