\name{extendTSS}
\alias{extendTSS}
\title{
Extend TSS
}
\description{
Extend TSS
}
\usage{
extendTSS(gene, seqlengths = NULL, genome = NULL,
    gene_id_type = NULL, mode = "basalPlusExt", extend_from = c("TSS", "gene"),
    basal_upstream = 5000, basal_downstream = 1000, extension = 1000000,
    verbose = great_opt$verbose, .attr = list())
}
\arguments{

  \item{gene}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object of gene (or TSS) coordinates.}
  \item{extend_from}{Should the gene be extended only from its TSS or the complete gene?}
  \item{seqlengths}{A named vector of chromosome lengths. If it is not provided, it is taken by \code{seqlengths(gene)}.}
  \item{genome}{UCSC genome can be set here, then \code{seqlengths} will be automatically retrieved from UCSC server.}
  \item{gene_id_type}{Gene ID types in \code{gene}. You need to set this argument if you use built-in gene sets in \code{\link{great}} so that genes can be correctly mapped. The value can only be one of "SYMBOL", "ENTREZ", "ENSEMBL" and "REFSEQ".}
  \item{mode}{The mode to extend TSS. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See "Details" section.}
  \item{basal_upstream}{In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.}
  \item{basal_downstream}{In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.}
  \item{extension}{Extensions from the basal domains. The value can also be a vector of length two which corresponds to extension to upstream and downstream respectively.}
  \item{verbose}{Whether to print messages.}
  \item{.attr}{Only used internally.}

}
\details{
Following are general explanations of the three modes for extending TSS:

\describe{
  \item{\code{basalPlusExt}}{1. TSS are extended into basal domains (e.g. by upstream 5kb, downstream 1kb); 2. basal domains are sorted by their genomic coordinates; 3. each basal domain is extended to its both sides until it reaches the next TSS's basal domain or it reaches the maximal extension (e.g. 1000kb).}
  \item{\code{twoClosest}}{1. TSS are sorted by their genomic coordinates; 2. each TSS is extended to its both sides until it reaches the next TSS or it reaches the maximal extension (e.g. 1000kb).}
  \item{\code{oneClosest}}{1. TSS are sorted by their genomic coordinates; 2. each TSS is extended to its both sides until it reaches the middle point of itself and the next TSS or it reaches the maximal extension (e.g. 1000kb).}
}

The official explanation is at \url{https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules} .
}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with one meta column 'gene_id'.
}
\examples{
# There is no example
NULL

}
