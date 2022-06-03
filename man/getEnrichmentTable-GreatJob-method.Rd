\name{getEnrichmentTable-GreatJob-method}
\alias{getEnrichmentTable,GreatJob-method}
\title{
Get a single enrichment table from GREAT web server
}
\description{
Get a single enrichment table from GREAT web server
}
\usage{
\S4method{getEnrichmentTable}{GreatJob}(object, ontology, ...)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}
  \item{ontology}{A single ontology names. Valid values are in \code{\link{availableOntologies}}. }
  \item{...}{All pass to \code{\link{getEnrichmentTables,GreatJob-method}}.}

}
\value{
A data frame of the enrichment results for a single ontology.
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "GreatJob.rds", package = "rGREAT"))
tb = getEnrichmentTable(job, ontology = "GO Molecular Function")
head(tb)
}
