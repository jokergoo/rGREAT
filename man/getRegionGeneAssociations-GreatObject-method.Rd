\name{getRegionGeneAssociations-GreatObject-method}
\alias{getRegionGeneAssociations,GreatObject-method}
\title{
Get region-gene associations
}
\description{
Get region-gene associations
}
\usage{
\S4method{getRegionGeneAssociations}{GreatObject}(object, term_id = NULL, use_symbols = TRUE)
}
\arguments{

  \item{object}{A \code{\link{GreatObject-class}} object returned by \code{\link{great}}.}
  \item{term_id}{Term ID.}
  \item{use_symbols}{Whether to use gene symbols}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. Please the two meta columns are in formats of \code{CharacterList}
and \code{IntegerList} because a region may associate to multiple genes.
}
\examples{
obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
getRegionGeneAssociations(obj)
}
