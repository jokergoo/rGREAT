\name{availableOntologies-GreatJob-method}
\alias{availableOntologies,GreatJob-method}
\alias{availableOntologies}
\title{
All available ontology names of the GREAT job
}
\description{
All available ontology names of the GREAT job
}
\usage{
\S4method{availableOntologies}{GreatJob}(object, category = NULL)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}
  \item{category}{one or multiple categories. All available categories can be got by \code{\link{availableCategories}}.}

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
