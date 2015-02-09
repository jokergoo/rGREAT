\name{GreatJob-class}
\docType{class}
\alias{GreatJob-class}
\alias{GreatJob}
\title{
Class to store and retrieve GREAT results  


}
\description{
Class to store and retrieve GREAT results  


}
\details{
After submitting request to GREAT server, the generated results will be  available on GREAT server for some time. The \code{GreatJob} class is defined  to store parameters that user has set and result tables what were retrieved  from GREAT server.  


}
\section{Constructor}{
Users don't need to construct by hand, \code{\link{submitGreatJob}} is used to generate a \code{GreatJob} instance.  


}
\section{Workflow}{
After submitting request to GREAT server, users can perform following steps:  

\itemize{
  \item call \code{\link{getEnrichmentTables}} to get enrichment tables for selected ontologies catalogues.
  \item call \code{\link{plotRegionGeneAssociationGraphs}} to get associations between regions and genes as well as making plots.  
}


}
\section{Constructor}{
Users don't need to construct by hand, \code{\link{submitGreatJob}} is used to generate a \code{GreatJob} instance.
}
\section{Workflow}{
After submitting request to GREAT server, users can perform following steps:  

\itemize{
  \item call \code{\link{getEnrichmentTables}} to get enrichment tables for selected ontologies catalogues.
  \item call \code{\link{plotRegionGeneAssociationGraphs}} to get associations between regions and genes as well as making plots.  
}
}
\author{
Zuguang gu <z.gu@dkfz.de>  


}
\examples{
# please refer to page of `submitGreatJob`
}
