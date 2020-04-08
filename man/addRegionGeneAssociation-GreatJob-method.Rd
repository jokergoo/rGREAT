\name{addRegionGeneAssociation-GreatJob-method}
\alias{addRegionGeneAssociation,GreatJob-method}
\title{
Add region-gene association to the internal enrichment tables
}
\description{
Add region-gene association to the internal enrichment tables
}
\usage{
\S4method{addRegionGeneAssociation}{GreatJob}(job, ontology = names(job@enrichment_tables), verbose = TRUE)
}
\arguments{

  \item{job}{a \code{\link{GreatJob-class}} instance}
  \item{ontology}{ontology name}
  \item{verbose}{whether show message}

}
\details{
After successfully run this function, users need to rerun \code{\link{getEnrichmentTables}} to 
get the enrichment tables that contains the new gene association column.
}
\alias{addRegionGeneAssociation}
\examples{
# There is no example
NULL
}
