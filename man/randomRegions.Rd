\name{randomRegions}
\alias{randomRegions}
\title{
Generate random regions
}
\description{
Generate random regions
}
\usage{
randomRegions(genome = NULL, nr = 1000, seqlengths = NULL,
    width_fun = function(n) runif(n, min = 1000, max = 10000))
}
\arguments{

  \item{genome}{UCSC genome version, e.g. "hg19".}
  \item{nr}{Number of regions.}
  \item{seqlengths}{Alternatively, you can also specify a named vector of seqlengths (chromosome lengths).}
  \item{width_fun}{A function which defines the distribution of region widths.}

}
\details{
The number of regions per chromosome is proportional to the chromsome length.
}
\examples{
gr = randomRegions(genome = "hg19")
quantile(width(gr))
}
