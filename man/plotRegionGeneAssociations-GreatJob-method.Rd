\name{plotRegionGeneAssociations-GreatJob-method}
\alias{plotRegionGeneAssociations,GreatJob-method}
\title{
Plot region-gene associations
}
\description{
Plot region-gene associations
}
\usage{
\S4method{plotRegionGeneAssociations}{GreatJob}(object, ontology = NULL, term_id = NULL, which_plot = 1:3,
    request_interval = 10, max_tries = 100, verbose = great_opt$verbose)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}
  \item{ontology}{A single ontology names. Valid values are in \code{\link{availableOntologies}}. }
  \item{term_id}{Term id in the selected ontology}
  \item{which_plot}{Which plots to draw? The value should be in \code{1, 2, 3}. See "Details" section for explanation.}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximal times for automatically reconnecting GREAT web server.}
  \item{verbose}{Whether to show messages.}

}
\details{
There are following figures:

\itemize{
  \item Association between regions and genes (\code{which_plot = 1}).
  \item Distribution of distance to TSS (\code{which_plot = 2}).
  \item Distribution of absolute distance to TSS (\code{which_plot = 3}).
}

If \code{ontology} and \code{term_id} are set, only regions and genes corresponding to 
selected ontology term will be used. Valid value for \code{ontology} is in 
\code{\link{availableOntologies}} and valid value for \code{term_id} is from 'id' column 
in the table which is returned by \code{\link{getEnrichmentTables}}.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))

plotRegionGeneAssociations(job)
plotRegionGeneAssociations(job, which_plot = 1)
plotRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")
}
