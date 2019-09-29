\name{availableCategories-GreatJob-method}
\alias{availableCategories,GreatJob-method}
\alias{availableCategories}
\title{
Available ontology categories
}
\description{
Available ontology categories
}
\usage{
\S4method{availableCategories}{GreatJob}(job)
}
\arguments{

  \item{job}{a \code{\link{GreatJob-class}} instance}

}
\details{
The values of the supported categories sometime change. You should run the function to get the real-time
values. The meaning of categories returned is quite self-explained by the name.
}
\value{
The returned value is a vector of categories.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
# note the `job` was generated from GREAT 3.0.0
job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
availableCategories(job)
}
