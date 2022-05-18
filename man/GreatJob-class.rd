\name{GreatJob-class}
\docType{class}
\alias{GreatJob-class}
\title{
Class to store and retrieve GREAT results
}
\description{
Class to store and retrieve GREAT results
}
\details{
After submitting request to GREAT server, the generated results will be 
available on GREAT server for some time. The \code{GreatJob-class} is defined 
to store parameters that user has set and result tables what were retrieved 
from GREAT server.
}
\section{Constructor}{
Users don't need to construct by hand, \code{\link{submitGreatJob}} is used to generate a \code{GreatJob-class} instance.}
\section{Workflow}{
After submitting request to GREAT server, users can perform following steps:

\itemize{
  \item \code{\link{getEnrichmentTables,GreatJob-method}} to get enrichment tables for selected ontologies catalogues.
  \item \code{\link{plotRegionGeneAssociations,GreatJob-method}} to plot associations between regions and genes
  \item \code{\link{getRegionGeneAssociations,GreatJob-method}} to get a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object which contains associations bewteen regions and genes.
  \item \code{\link{shinyReport,GreatJob-method}} to view the results by a shiny application.
}}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
