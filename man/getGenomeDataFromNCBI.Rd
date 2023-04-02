\name{getGenomeDataFromNCBI}
\alias{getGenomeDataFromNCBI}
\title{
Get genome data from NCBI
}
\description{
Get genome data from NCBI
}
\usage{
getGenomeDataFromNCBI(refseq_assembly_accession, return_granges = FALSE)
}
\arguments{

  \item{refseq_assembly_accession}{The RefSeq accession number for the assembly, such as "GCF_000001405.40" for human.}
  \item{return_granges}{If the assembly is already on chromosome level, it will directly construct a GRanges object where "chromosomes" are only used and chromosome lengths are corrected fitted in its \code{seqlengths}. }

}
\details{
Only protein coding genes are used.
}
\value{
If \code{return_granges} is set to \code{FALSE}, it returns a list of two data frames:

\describe{
  \item{genome}{A data frame of several columns.}
  \item{gene}{A data frame for genes. The first column contains the RefSeq accession numbers of the corresponding contigs. If the genome is assembled on the chromosome level, the first column corresponds to chromosomes. The contig names can be converted to other names with the information in the \code{genome} data frame.}
}
}
\examples{
if(FALSE) {
getGenomeDataFromNCBI("GCF_000001405.40", return_granges = TRUE)
getGenomeDataFromNCBI("GCF_000001405.40")
}
}
