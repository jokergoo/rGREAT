\name{reduce_by_start_and_end}
\alias{reduce_by_start_and_end}
\title{
Reduce by start and end
}
\description{
Reduce by start and end
}
\usage{
reduce_by_start_and_end(s, e)
}
\arguments{

  \item{s}{Start positions. Sorted.}
  \item{e}{End positions. Sorted.}

}
\details{
Only internally used.
}
\value{
Sum of total widths of the reduced regions.
}
\examples{
if(FALSE) {
getGenomeDateFromNCBI("GCF_000001405.40", return_granges = TRUE)
getGenomeDateFromNCBI("GCF_000001405.40")
}

}
