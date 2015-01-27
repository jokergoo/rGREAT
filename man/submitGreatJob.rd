\name{submitGreatJob}
\alias{submitGreatJob}
\title{
  Send requests to GREAT web server  


}
\description{
  Send requests to GREAT web server  


}
\usage{
submitGREATJob(gr, bg = NULL,
    species               = c("hg19", "hg18", "mm9", "danRer7"),
    includeCuratedRegDoms = TRUE,
    bgChoice              = c("wholeGenome", "data"),
    rule                  = c("basalPlusExt", "twoClosest", "oneClosest"),
    adv_upstream          = 5.0,
    adv_downstream        = 1.0,
    adv_span              = 1000.0,
    adv_twoDistance       = 1000.0,
    adv_oneDistance       = 1000.0,
    request_interval = 300,
    max_tries = 10)
}
\arguments{
  \item{gr}{A \code{\link[GenomicRanges]{GRanges}} object or a data frame which contains at least three columns (chr, start and end). Regions for test.}
  \item{bg}{A \code{\link[GenomicRanges]{GRanges}} object or a data frame. Background regions if needed.}
  \item{species}{Species. Only four species ("hg19", "hg18", "mm9", "danRer7") are supported.}
  \item{includeCuratedRegDoms}{Whether to include curated regulatory domains.}
  \item{bgChoice}{How to define background. If it is set as \code{data}, \code{bg} should be set as well.}
  \item{rule}{How to associate genomic regions to genes. See 'details' section.}
  \item{adv_upstream}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_downstream}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_span}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_twoDistance}{Unit: kb, only used when rule is \code{twoClosest}}
  \item{adv_oneDistance}{Unit: kb, only used when rule is \code{oneClosest}}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximum times trying to connect to GREAT web server.}

}
\details{
  Note it is not the standard GREAT API. This function directly send data to GREAT web server by HTTP POST.  

  Following text is copied from GREAT web site ( \url{http://bejerano-test.stanford.edu/great/public/html/index.php} )  

  Explanation of \code{rule} and settings with names started with 'adv_':  

\describe{
  \item{basalPlusExt}{Mode 'Basal plus extension'. Gene regulatory domain definition:  Each gene is assigned a basal regulatory domain of a minimum distance upstream  and downstream of the TSS (regardless of other nearby genes, controlled by \code{adv_upstream} and  \code{adv_downstream} argument). The gene regulatory domain is extended in both directions  to the nearest gene's basal domain but no more than the maximum extension in one direction (controlled by \code{adv_span}).}
  \item{twoClosest}{Mode 'Two nearest genes'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the nearest  gene's TSS (controlled by \code{adv_twoDistance}) but no more than the maximum extension in one direction.}
  \item{oneClosest}{Mode 'Single nearest gene'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the midpoint  between the gene's TSS and the nearest gene's TSS (controlled by \code{adv_oneDistance}) but no more than the maximum  extension in one direction.}
}

}
\value{
  A \code{GREAT_Job} class object which can be used to get results from GREAT server.  


}
\seealso{
  \code{\link{GREAT_Job-class}}  


}
\author{
  Zuguang gu <z.gu@dkfz.de>  


}
\examples{
# see `GREAT_Job` man page or the vignette
}