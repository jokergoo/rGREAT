\name{plotRegionGeneAssociationGraphs}
\alias{plotRegionGeneAssociationGraphs}
\title{
  Plot region-gene association figures  


}
\description{
  Plot region-gene association figures  


}
\usage{
plotRegionGeneAssociationGraphs(job, type = 1:3, ontology = NULL,
    termID = NULL, request_interval = 30, max_tries = 100)
}
\arguments{
  \item{job}{job object which is returned by \code{\link{submitGreatJob}} }
  \item{type}{type of plots, should be in \code{1, 2, 3}}
  \item{ontology}{ontology name}
  \item{termID}{term id}
  \item{request_interval}{interval for two requests. Default is 300 seconds.}
  \item{max_tries}{maximum tries}

}
\details{
  Figures are:  

  \itemize{
    \item association between regions and genes
    \item distribution of distance to TSS
    \item distribution of absolute distance to TSS
  }
  If \code{ontology} and \code{termID} are set, only regions and genes corresponding to selected ontology term will be used. Valid value for \code{ontology} is in \code{\link{availableOntologies}} and valid value for \code{termID} is from 'id' column in the table which is returned by \code{\link{getGreatTable}}.  


}
\value{
  a \code{GRanges} object. Columns in metadata are:  

\describe{
  \item{gene}{genes that are associated with corresponding regions}
  \item{distTSS}{distance from the regions to TSS of the associated gene}
}
  The returned values corresponds to whole input regions or only regions in specified ontology term, depending on user's settings. 


}
\examples{
set.seed(123)
bed = circlize:::generateRandomBed(nr = 1000, nc = 0)
	
#job = submitGreatJob(bed)
load(paste0(system.file(package = "rGREAT"), "/extdata/job.RData"))

df = plotRegionGeneAssociationGraphs(job)
df = plotRegionGeneAssociationGraphs(job, type = 2)

df = plotRegionGeneAssociationGraphs(job, type = 1,
    onto = "GO_Molecular_Function", termID = "GO:0004984")
}
