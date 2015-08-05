\name{availableOntologies-GreatJob-method}
\alias{availableOntologies,GreatJob-method}
\alias{availableOntologies}
\title{
All available ontology names

}
\description{
All available ontology names

}
\usage{
\S4method{availableOntologies}{GreatJob}(job, category = NULL)}
\arguments{

  \item{job}{a \code{\link{GreatJob}} instance}
  \item{category}{one or multiple categories. All available categories can be get by \code{\link{availableCategories}}}
}
\details{
The values of the supported ontologies sometime change. You should run the function to get the realtime
values. The meaning of ontology returned is quite self-explained by the name.

}
\value{
The returned values is a vector of ontologies.

}
\author{
Zuguang gu <z.gu@dkfz.de>

}
\examples{

job = readRDS(paste0(system.file("extdata", package = "rGREAT"), "/job.rds"))
availableOntologies(job)
availableOntologies(job, category = "Pathway Data")
}
