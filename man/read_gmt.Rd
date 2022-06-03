\name{read_gmt}
\alias{read_gmt}
\title{
Read gmt gene sets file
}
\description{
Read gmt gene sets file
}
\usage{
read_gmt(x, from = NULL, to = NULL, orgdb = NULL)
}
\arguments{

  \item{x}{The file name of a .gmt file.}
  \item{from}{Gene ID type in the original gmt file. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.}
  \item{to}{Gene ID type that you want to convert to. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.}
  \item{orgdb}{The name of an OrgDb database.}

}
\value{
A named list of vectors.
}
\examples{
read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"))
}
