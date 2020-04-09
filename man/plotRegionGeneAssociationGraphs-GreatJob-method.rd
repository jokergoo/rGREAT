\name{plotRegionGeneAssociationGraphs-GreatJob-method}
\alias{plotRegionGeneAssociationGraphs,GreatJob-method}
\alias{plotRegionGeneAssociationGraphs}
\title{
Plot region-gene association figures
}
\description{
Plot region-gene association figures
}
\usage{
\S4method{plotRegionGeneAssociationGraphs}{GreatJob}(job, type = 1:3, ontology = NULL,
    termID = NULL, request_interval = 30, max_tries = 100, verbose = TRUE,
    plot = TRUE)
}
\arguments{

  \item{job}{a \code{\link{GreatJob-class}} instance}
  \item{type}{type of plots, should be in \code{1, 2, 3}. See details section for explanation}
  \item{ontology}{ontology name}
  \item{termID}{term id which corresponds to the selected ontology}
  \item{request_interval}{time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{maximum tries}
  \item{verbose}{whether show message}
  \item{plot}{whether make plots}

}
\details{
Generated figures are:

\itemize{
  \item association between regions and genes
  \item distribution of distance to TSS
  \item distribution of absolute distance to TSS
}

If \code{ontology} and \code{termID} are set, only regions and genes corresponding to 
selected ontology term will be used. Valid value for \code{ontology} is in 
\code{\link{availableOntologies}} and valid value for \code{termID} is from 'id' column 
in the table which is returned by \code{\link{getEnrichmentTables}}.
}
\value{
a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. Columns in metadata are:

\describe{
  \item{gene}{genes that are associated with corresponding regions}
  \item{distTSS}{distance from the regions to TSS of the associated gene}
}

The returned values corresponds to whole input regions or only regions in specified ontology term, 
depending on user's setting.

If there is no gene associated with the region, corresponding \code{gene} and \code{distTSS}
columns will be \code{NA}.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))

res = plotRegionGeneAssociationGraphs(job)
res

plotRegionGeneAssociationGraphs(job, type = 1)

res = plotRegionGeneAssociationGraphs(job, ontology = "GO Molecular Function",
    termID = "GO:0004984")
res
}
