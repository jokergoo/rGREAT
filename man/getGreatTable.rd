\name{getGreatTable}
\alias{getGreatTable}
\title{
  Get enrichment tables from GREAT web server  


}
\description{
  Get enrichment tables from GREAT web server  


}
\usage{
getGreatTable(job, ontology = NULL, category = c("GO", "Pathway_Data"),
    request_interval = 30, max_tries = 100)
}
\arguments{
  \item{job}{job object which is returned by \code{\link{submitGreatJob}}}
  \item{ontology}{ontology names. Valid values are in \code{\link{availableOntologies}}. \code{ontology} is prior to \code{category} argument.}
  \item{category}{Pre-defined categories. Valid values are in \code{\link{availableCategories}}}
  \item{request_interval}{interval for two requests. Default is 300 seconds.}
  \item{max_tries}{maximum tries}

}
\details{
  Please note there is no FDR column in original tables. Users should calculate by themselves by functions such as \code{\link[stats]{p.adjust}}  


}
\value{
  a list of data frames which are as same as in GREAT website. 


}
\examples{
set.seed(123)
bed = circlize:::generateRandomBed(nr = 1000, nc = 0)
	
#job = submitGreatJob(bed)
load(paste0(system.file(package = "rGREAT"), "/extdata/job.RData"))

tb = getGreatTable(job)
names(tb)
head(tb[[1]])
}
