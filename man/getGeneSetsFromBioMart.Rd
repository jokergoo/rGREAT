\name{getGeneSetsFromBioMart}
\alias{getGeneSetsFromBioMart}
\title{
Get gene sets from BioMart
}
\description{
Get gene sets from BioMart
}
\usage{
getGeneSetsFromBioMart(dataset, ontology = "bp")
}
\arguments{

  \item{dataset}{Name of the dataset.}
  \item{ontology}{Value should be bp, mf or cc.}

}
\details{
GO gene sets are from \code{\link[BioMartGOGeneSets]{getBioMartGOGeneSets}}.
}
\value{
A list of vectors where each vector contains Ensembl IDs annotated to a GO term.
}
\examples{
# There is no example
NULL

}
