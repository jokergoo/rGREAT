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
    bgChoice              = ifelse(is.null(bg), "wholeGenome", "data"),
    rule                  = c("basalPlusExt", "twoClosest", "oneClosest"),
    adv_upstream          = 5.0,
    adv_downstream        = 1.0,
    adv_span              = 1000.0,
    adv_twoDistance       = 1000.0,
    adv_oneDistance       = 1000.0,
    request_interval = 300,
    max_tries = 10,
    version = "default",
    base_url = "http://great.stanford.edu/public/cgi-bin/")
}
\arguments{

  \item{gr}{A \code{\link[GenomicRanges]{GRanges}} object or a data frame which contains at least three columns (chr, start and end). Regions for test.}
  \item{bg}{A \code{\link[GenomicRanges]{GRanges}} object or a data frame. Background regions if needed.}
  \item{species}{Species. "hg19", "mm10", "mm9", "danRer7" are supported in GREAT version 3.x.x and "hg19", "hg18", "mm9", "danRer7" are supported in GREAT version 2.x.x.}
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
  \item{version}{version of GREAT. The value should be "3.0.0", "2.0.2". Shorten version numbers can also be used, such as using "3" or "3.0" is same as "3.0.0".}
  \item{base_url}{the url of \code{cgi-bin} path, only used when explicitly specified.}

}
\details{
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

When \code{bg} is set, some pre-processing is applied before submitting to GREAT server for the reason
that GREAT needs \code{gr} should be exactly subsets of \code{bg}, which means for any region in \code{gr}, there
must be a region in \code{bg} which is exactly the same. Taking following example:

for \code{gr}:

  \preformatted{
    chr1 200 300
    chr1 250 400  }

for \code{bg}:

  \preformatted{
    chr1 100 250
    chr1 300 500
    chr1 400 600  }

They will be transformed as: for \code{gr}:

  \preformatted{
    chr1 200 250
    chr1 300 400  }

for \code{bg}:

  \preformatted{
    chr1 100 199
    chr1 200 250
    chr1 300 400
    chr1 401 600  }
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
job = submitGreatJob(bed)

# more parameters can be set for the job
\dontrun{
job = submitGreatJob(bed, species = "mm9")
job = submitGreatJob(bed, bg, species = "mm9", bgChoise = "data")
job = submitGreatJob(bed, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
job = submitGreatJob(bed, rule = "twoClosest", adv_twoDistance = 2000)
job = submitGreatJob(bed, rule = "oneClosest", adv_oneDistance = 2000)
}

}
