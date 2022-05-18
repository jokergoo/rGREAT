\name{plotRegionGeneAssociations-GreatObject-method}
\alias{plotRegionGeneAssociations,GreatObject-method}
\title{
Plot region-gene associations
}
\description{
Plot region-gene associations
}
\usage{
\S4method{plotRegionGeneAssociations}{GreatObject}(object, term_id = NULL, which_plot = 1:3)
}
\arguments{

  \item{object}{A \code{\link{GreatObject-class}} object returned by \code{\link{great}}.}
  \item{term_id}{Term ID.}
  \item{which_plot}{Which plots to draw? The value should be in \code{1, 2, 3}. See "Details" section for explanation.}

}
\details{
There are following figures:

\itemize{
  \item Association between regions and genes (\code{which_plot = 1}).
  \item Distribution of distance to TSS (\code{which_plot = 2}).
  \item Distribution of absolute distance to TSS (\code{which_plot = 3}).
}
}
\examples{
obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
plotRegionGeneAssociations(obj)
}
