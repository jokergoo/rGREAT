\name{getRegionGeneAssociations-GreatJob-method}
\alias{getRegionGeneAssociations,GreatJob-method}
\title{
Get region-gene associations
}
\description{
Get region-gene associations
}
\usage{
\S4method{getRegionGeneAssociations}{GreatJob}(object, ontology = NULL, term_id = NULL,
    request_interval = 10, max_tries = 100, verbose = great_opt$verbose)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}
  \item{ontology}{ontology name}
  \item{term_id}{Term id in the selected ontology.}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximal times for automatically reconnecting GREAT web server.}
  \item{verbose}{Whether to show messages.}

}
\value{
A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. Please the two meta columns are in formats of \code{CharacterList}
and \code{IntegerList} because a region may associate to multiple genes.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))

gr = getRegionGeneAssociations(job)
gr
}
