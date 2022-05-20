\name{getEnrichmentTable-GreatObject-method}
\alias{getEnrichmentTable,GreatObject-method}
\title{
Get enrichment table
}
\description{
Get enrichment table
}
\usage{
\S4method{getEnrichmentTable}{GreatObject}(object, min_region_hits = 5)
}
\arguments{

  \item{object}{A \code{\link{GreatObject-class}} object returned by \code{\link{great}}.}
  \item{min_region_hits}{Minimal number of input regions overlapping to the geneset associated regions.}

}
\details{
Note: adjusted p-values are re-calculated based on \code{min_region_hits}.
}
\value{
A data frame of enrichment results
}
\examples{
obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
getEnrichmentTable(obj)
}
