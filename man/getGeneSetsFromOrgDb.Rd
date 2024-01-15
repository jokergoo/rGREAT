\name{getGeneSetsFromOrgDb}
\alias{getGeneSetsFromOrgDb}
\title{
Get GO gene sets from OrgDb object
}
\description{
Get GO gene sets from OrgDb object
}
\usage{
getGeneSetsFromOrgDb(orgdb, ontology = "BP")
}
\arguments{

  \item{orgdb}{An \code{OrgDb} object.}
  \item{ontology}{Value should be bp, mf or cc.}

}
\value{
A list of vectors where each vector contains Entrez IDs annotated to a GO term.
}
\examples{
# There is no example
NULL

}
