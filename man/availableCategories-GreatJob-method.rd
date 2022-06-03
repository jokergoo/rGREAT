\name{availableCategories-GreatJob-method}
\alias{availableCategories,GreatJob-method}
\alias{availableCategories}
\title{
Available ontology categories of the GREAT job
}
\description{
Available ontology categories of the GREAT job
}
\usage{
\S4method{availableCategories}{GreatJob}(object)
}
\arguments{

  \item{object}{A \code{\link{GreatJob-class}} object returned by \code{\link{submitGreatJob}}.}

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
job = readRDS(system.file("extdata", "GreatJob.rds", package = "rGREAT"))
availableCategories(job)
}
