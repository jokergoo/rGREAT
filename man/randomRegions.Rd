\name{randomRegions}
\alias{randomRegions}
\title{
Generate random regions
}
\description{
Generate random regions
}
\usage{
randomRegions(nr = 1000, genome = NULL,
    width_fun = function(n) runif(n, min = 1000, max = 10000))
}
\arguments{

  \item{nr}{Number of regions.}
  \item{genome}{UCSC genome version, e.g. "hg19".}
  \item{width_fun}{A function which defines the distribution of region widths.}

}
\details{
The number of regions per chromosome is proportional to the chromsome length.
}
\examples{
gr = randomRegions(1000)
quantile(width(gr))
}
