\name{extendTSSFromOrgDb}
\alias{extendTSSFromOrgDb}
\title{
Extend TSS
}
\description{
Extend TSS
}
\usage{
extendTSSFromOrgDb(orgdb, verbose = great_opt$verbose, ...)
}
\arguments{

  \item{orgdb}{Name of "org.*" packages from Bioconductor. All supported OrgDb packages are in \code{rGREAT:::BIOC_ANNO_PKGS$orgdb}. }
  \item{verbose}{Whether to print messages.}
  \item{...}{All pass to \code{\link{extendTSS}}.}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with one meta column 'gene_id'.
}
\examples{
if(FALSE) {
extendTSSFromOrgDb("Org.Hs.eg.db")
extendTSSFromOrgDb("hg19")
}
}
