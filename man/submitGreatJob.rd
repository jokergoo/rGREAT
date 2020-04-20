\name{submitGreatJob}
\alias{submitGreatJob}
\title{
Send requests to GREAT web server
}
\description{
Send requests to GREAT web server
}
\usage{
submitGreatJob(gr, bg = NULL,
    species               = "hg19",
    includeCuratedRegDoms = TRUE,
    rule                  = c("basalPlusExt", "twoClosest", "oneClosest"),
    adv_upstream          = 5.0,
    adv_downstream        = 1.0,
    adv_span              = 1000.0,
    adv_twoDistance       = 1000.0,
    adv_oneDistance       = 1000.0,
    request_interval = 60,
    max_tries = 10,
    version = DEFAULT_VERSION,
    base_url = "http://great.stanford.edu/public/cgi-bin")
}
\arguments{

  \item{gr}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object or a data frame which contains at least three columns (chr, start and end). Regions for test.}
  \item{bg}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object or a data frame. Background regions if needed. Note \code{gr} should be exactly subset of \code{bg} for all columns in \code{gr}. Check \url{http://great.stanford.edu/help/display/GREAT/File+Formats#FileFormats-Whatshouldmybackgroundregionsfilecontain\%3F} for more explanation.}
  \item{species}{Species. "hg38", "hg19", "mm10", "mm9" are supported in GREAT version 4.x.x, "hg19", "mm10", "mm9", "danRer7" are supported in GREAT version 3.x.x and "hg19", "hg18", "mm9", "danRer7" are supported in GREAT version 2.x.x.}
  \item{includeCuratedRegDoms}{Whether to include curated regulatory domains.}
  \item{rule}{How to associate genomic regions to genes. See 'details' section.}
  \item{adv_upstream}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_downstream}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_span}{Unit: kb, only used when rule is \code{basalPlusExt}}
  \item{adv_twoDistance}{Unit: kb, only used when rule is \code{twoClosest}}
  \item{adv_oneDistance}{Unit: kb, only used when rule is \code{oneClosest}}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximum times trying to connect to GREAT web server.}
  \item{version}{version of GREAT. The value should be "4.0.4", "3.0.0", "2.0.2". Shorten version numbers can also be used, such as using "4" or "4.0" is same as "4.0.4".}
  \item{base_url}{the url of \code{cgi-bin} path, only used when explicitly specified.}

}
\details{
Note: [On Aug 19 2019 GREAT released version 4](\url{http://great.stanford.edu/help/display/GREAT/Version+History}  where it supports \code{hg38} genome and removes some ontologies such pathways. \code{\link{submitGreatJob}} still
takes \code{hg19} as default. \code{hg38} can be specified by the \code{species = "hg38"} argument.
To use the older versions such as 3.0.0, specify as \code{submitGreatJob(..., version = "3.0.0")}.

Note it is not the standard GREAT API. This function directly send data to GREAT web server
by HTTP POST.

Following text is copied from GREAT web site ( \url{http://great.stanford.edu/public/html/} )

Explanation of \code{rule} and settings with names started with 'adv_' (advanced settings):

\describe{
  \item{basalPlusExt}{Mode 'Basal plus extension'. Gene regulatory domain definition:  Each gene is assigned a basal regulatory domain of a minimum distance upstream  and downstream of the TSS (regardless of other nearby genes, controlled by \code{adv_upstream} and  \code{adv_downstream} argument). The gene regulatory domain is extended in both directions  to the nearest gene's basal domain but no more than the maximum extension in one direction (controlled by \code{adv_span}).}
  \item{twoClosest}{Mode 'Two nearest genes'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the nearest  gene's TSS (controlled by \code{adv_twoDistance}) but no more than the maximum extension in one direction.}
  \item{oneClosest}{Mode 'Single nearest gene'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the midpoint  between the gene's TSS and the nearest gene's TSS (controlled by \code{adv_oneDistance}) but no more than the maximum  extension in one direction.}
}
}
\section{Note}{
takes \code{hg19} as default. \code{hg38} can be specified by the \code{species = "hg38"} argument.
To use the older versions such as 3.0.0, specify as \code{submitGreatJob(..., version = "3.0.0")}.**

Note it is not the standard GREAT API. This function directly send data to GREAT web server
by HTTP POST.

Following text is copied from GREAT web site ( \url{http://great.stanford.edu/public/html/} )

Explanation of \code{rule} and settings with names started with 'adv_' (advanced settings):

\describe{
  \item{basalPlusExt}{Mode 'Basal plus extension'. Gene regulatory domain definition:  Each gene is assigned a basal regulatory domain of a minimum distance upstream  and downstream of the TSS (regardless of other nearby genes, controlled by \code{adv_upstream} and  \code{adv_downstream} argument). The gene regulatory domain is extended in both directions  to the nearest gene's basal domain but no more than the maximum extension in one direction (controlled by \code{adv_span}).}
  \item{twoClosest}{Mode 'Two nearest genes'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the nearest  gene's TSS (controlled by \code{adv_twoDistance}) but no more than the maximum extension in one direction.}
  \item{oneClosest}{Mode 'Single nearest gene'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the midpoint  between the gene's TSS and the nearest gene's TSS (controlled by \code{adv_oneDistance}) but no more than the maximum  extension in one direction.}
}
}
\value{
A \code{\link{GreatJob-class}} class object which can be used to get results from GREAT server.
}
\seealso{
\code{\link{GreatJob-class}}
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
bed = circlize::generateRandomBed(nr = 1000, nc = 0)
job = submitGreatJob(bed, version = "3.0.0")
job

# more parameters can be set for the job
if(FALSE) { # suppress running it when building the package
    # current GREAT version is 4.0.1
    job = submitGreatJob(bed, species = "mm9")
    job = submitGreatJob(bed, bg, species = "mm9", bgChoise = "data")
    job = submitGreatJob(bed, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
    job = submitGreatJob(bed, rule = "twoClosest", adv_twoDistance = 2000)
    job = submitGreatJob(bed, rule = "oneClosest", adv_oneDistance = 2000)
}
}
