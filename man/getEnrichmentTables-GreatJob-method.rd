\name{getEnrichmentTables-GreatJob-method}
\alias{getEnrichmentTables,GreatJob-method}
\alias{getEnrichmentTables}
\title{
Get enrichment tables from GREAT web server
}
\description{
Get enrichment tables from GREAT web server
}
\usage{
\S4method{getEnrichmentTables}{GreatJob}(job, ontology = NULL, category = "GO",
    request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
    verbose = TRUE)
}
\arguments{

  \item{job}{a \code{\link{GreatJob-class}} instance}
  \item{ontology}{ontology names. Valid values are in \code{\link{availableOntologies}}. \code{ontology} is prior to  \code{category} argument.}
  \item{category}{Pre-defined ontology categories. One category can contain more than one ontologies. Valid values are in  \code{\link{availableCategories}}}
  \item{request_interval}{time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{maximum tries}
  \item{download_by}{Internally used.}
  \item{verbose}{Whether print messages.}

}
\details{
The table contains statistics for the each term in each ontology catalogue.

Please note there is no FDR column in original tables. Users should 
calculate by themselves by functions such as \code{\link[stats]{p.adjust}}
}
\value{
The returned value is a list of data frames in which each one corresponds to 
result for a single ontology. The structure of the data frames are same as 
the tables available on GREAT website.
}
\section{See}{
\code{\link{availableOntologies}}, \code{\link{availableCategories}}}
\seealso{
\code{\link{availableOntologies}}, \code{\link{availableCategories}}

}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
tb = getEnrichmentTables(job)
names(tb)
head(tb[[1]])
job

tb = getEnrichmentTables(job, ontology = "GO Molecular Function")
tb = getEnrichmentTables(job, category = "GO")
}
