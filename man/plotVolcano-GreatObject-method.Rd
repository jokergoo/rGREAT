\name{plotVolcano-GreatObject-method}
\alias{plotVolcano,GreatObject-method}
\title{
Make volcano plot
}
\description{
Make volcano plot
}
\usage{
\S4method{plotVolcano}{GreatObject}(object, min_region_hits = 5,
    x_values = c("fold_enrichment", "z-score"),
    y_values = c("p_value", "p_adjust"),
    main = NULL)
}
\arguments{

  \item{object}{A \code{\link{GreatObject-class}} object returned by \code{\link{great}}.}
  \item{min_region_hits}{Minimal number of input regions overlapping to the geneset associated regions.}
  \item{x_values}{Which values for the x-axis.}
  \item{y_values}{Which values for the y-axis.}
  \item{main}{Title of the plot.}

}
\details{
Since the enrichment is an over-representation test, it is only the half volcano.
}
\examples{
# There is no example
NULL

}
