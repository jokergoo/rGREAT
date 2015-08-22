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
The values of the supported categories sometime change. You should run the function to get the realtime
values. The meaning of categories returned is quite self-explained by the name.
}
\value{
The returned value is a vector of categories.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
job = readRDS(paste0(system.file("extdata", package = "rGREAT"), "/job.rds"))
availableCategories(job)

}
