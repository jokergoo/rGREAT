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
\S4method{availableOntologies}{GreatJob}(job, category = NULL)
}
\arguments{

  \item{job}{a \code{\link{GreatJob-class}} instance}
  \item{category}{one or multiple categories. All available categories can be get by \code{\link{availableCategories}}}

}
\details{
The values of the supported ontologies sometime change. You should run the function to get the real-time
values. The meaning of ontology returned is quite self-explained by the name.
}
\value{
The returned values is a vector of ontologies.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
availableOntologies(job)
availableOntologies(job, category = "Pathway Data")
}
