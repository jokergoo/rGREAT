\name{getEnrichmentTables-GreatJob-method}
\alias{getEnrichmentTables,GreatJob-method}
\title{
Get enrichment tables from GREAT web server
}
\description{
Get enrichment tables from GREAT web server
}
\usage{
\S4method{getEnrichmentTables}{GreatJob}(object, ontology = NULL, category = "GO",
    request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
    verbose = TRUE)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}
  \item{ontology}{Ontology names. Valid values are in \code{\link{availableOntologies}}. \code{ontology} is prior to  \code{category} argument.}
  \item{category}{Pre-defined ontology categories. One category can contain more than one ontologies. Valid values are in  \code{\link{availableCategories}}}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximal times for automatically reconnecting GREAT web server.}
  \item{download_by}{Internally used.}
  \item{verbose}{Whether to print messages.}

}
\value{
The structure of the data frames are same as the tables available on GREAT website.
}
\section{See}{
\code{\link{availableOntologies}}, \code{\link{availableCategories}}}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "GreatJob.rds", package = "rGREAT"))
tbl = getEnrichmentTables(job)
names(tbl)
head(tbl[[1]])
job

tbl = getEnrichmentTables(job, ontology = "GO Molecular Function")
tbl = getEnrichmentTables(job, category = "GO")
}
