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
  \item{download_by}{Internally used. The complete enrichment table is provided as json data on the website, but there is no information of gene-region association. By setting \code{download_by = 'tsv'}, another URL from GREAT will be envoked which also contains detailed information of which genes are associated with each input region, but due to the size of the output, only top 500 terms will be returned. So if you do not really want the gene-region association column, take the default value of this argument. The columns that contain statistics are identical.}
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
job = readRDS(system.file("extdata", "GreatJob.rds", package = "rGREAT"))
tbl = getEnrichmentTables(job)
names(tbl)
head(tbl[[1]])
job

tbl = getEnrichmentTables(job, ontology = "GO Molecular Function")
tbl = getEnrichmentTables(job, category = "GO")
}
