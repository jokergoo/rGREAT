\name{submitGreatJob}
\alias{submitGreatJob}
\title{
  Automatically send requests to GREAT web server  


}
\description{
  Automatically send requests to GREAT web server  


}
\usage{
submitGreatJob(gr, bg = NULL,
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
  \item{gr}{GRanges object or a data frame which contains at least three columns (chr, start and end). Test regions.}
  \item{bg}{GRanges object or a data frame. Background regions if needed.}
  \item{species}{Species. Only four species are supported.}
  \item{includeCuratedRegDoms}{Whether to include curated regulatory domains.}
  \item{bgChoice}{How to define background. If it is set as \code{data}, \code{bg} should be set as well.}
  \item{rule}{How to associate genomic regions to genes. See details section.}
  \item{adv_upstream}{unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_downstream}{unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_span}{unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_twoDistance}{unit: kb, only used when rule is \code{twoClosest}}
  \item{adv_oneDistance}{unit: kb, only used when rule is \code{oneClosest}}
  \item{request_interval}{interval for two requests. Default is 300 seconds.}
  \item{max_tries}{maximum tries}

}
\details{
  Note it is not the standard GREAT API. This function directly send data to GREAT web server by HTTP POST.  

  Following are copied from GREAT web site ( \url{http://bejerano-test.stanford.edu/great/public/html/index.php} )  

  Explanation of \code{rule}  

\describe{
  \item{basalPlusExt}{mode 'Basal plus extension'. Gene regulatory domain definition: Each gene is assigned a basal regulatory domain of a minimum distance upstream and downstream of the TSS (regardless of other nearby genes). The gene regulatory domain is extended in both directions to the nearest gene's basal domain but no more than the maximum extension in one direction.}
  \item{twoClosest}{mode 'Two nearest genes'. Gene regulatory domain definition: Each gene is assigned a regulatory domain that extends in both directions to the nearest gene's TSS but no more than the maximum extension in one direction.}
  \item{oneClosest}{mode 'Single nearest gene'. Gene regulatory domain definition: Each gene is assigned a regulatory domain that extends in both directions to the midpoint between the gene's TSS and the nearest gene's TSS but no more than the maximum extension in one direction.}
}

}
\value{
  a \code{GREAT_Job} object 


}
\examples{
set.seed(123)
bed = circlize:::generateRandomBed(nr = 1000, nc= 0)
	
#job = submitGreatJob(bed)
load(paste0(system.file(package = "rGREAT"), "/extdata/job.RData"))

#job = submitGreatJob(bed, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
#job = submitGreatJob(bed, rule = "twoClosest", adv_twoDistance = 2000)
#job = submitGreatJob(bed, rule = "oneClosest", adv_oneDistance = 2000)
}
