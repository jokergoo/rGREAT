\name{shinyReport-GreatJob-method}
\alias{shinyReport,GreatJob-method}
\title{
Shiny app on the GreatJob object
}
\description{
Shiny app on the GreatJob object
}
\usage{
\S4method{shinyReport}{GreatJob}(object)
}
\arguments{

  \item{object}{The \code{GreatJob} object returned by \code{\link{submitGreatJob}}.}

}
\value{
A shiny app object.
}
\examples{
if(FALSE) {
# pseudo code
job = submitGreatJob(...)
shinyReport(job)
}
}
