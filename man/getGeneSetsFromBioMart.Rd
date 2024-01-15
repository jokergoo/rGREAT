\name{getGeneSetsFromBioMart}
\alias{getGeneSetsFromBioMart}
\title{
Get GO gene sets from BioMart
}
\description{
Get GO gene sets from BioMart
}
\usage{
getGeneSetsFromBioMart(dataset, ontology = "bp")
}
\arguments{

  \item{dataset}{Name of the dataset.}
  \item{ontology}{Value should be bp, mf or cc.}

}
\details{
GO gene sets are from \code{BioMartGOGeneSets::getBioMartGOGeneSets}.
}
\value{
A list of vectors where each vector contains Ensembl IDs annotated to a GO term.
}
\examples{
# There is no example
NULL

}
