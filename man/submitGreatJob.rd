\name{submitGreatJob}
\alias{submitGreatJob}
\title{
Perform online GREAT analysis
}
\description{
Perform online GREAT analysis
}
\usage{
submitGreatJob(gr, bg = NULL,
    gr_is_zero_based      = FALSE,
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
    base_url = "http://great.stanford.edu/public/cgi-bin",
    help = TRUE)
}
\arguments{

  \item{gr}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object or a data frame which contains at least three columns (chr, start and end).}
  \item{bg}{Not supported any more. See explanations in section "When_background_regions_are_set".}
  \item{gr_is_zero_based}{Are start positions in \code{gr} zero-based?}
  \item{species}{Species. "hg38", "hg19", "mm10", "mm9" are supported in GREAT version 4.x.x, "hg19", "mm10", "mm9", "danRer7" are supported in GREAT version 3.x.x and "hg19", "hg18", "mm9", "danRer7" are supported in GREAT version 2.x.x.}
  \item{includeCuratedRegDoms}{Whether to include curated regulatory domains, see \url{https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules#AssociationRules-CuratedRegulatoryDomains} .}
  \item{rule}{How to associate genomic regions to genes. See 'Details' section.}
  \item{adv_upstream}{Unit: kb, only used when rule is \code{basalPlusExt}.}
  \item{adv_downstream}{Unit: kb, only used when rule is \code{basalPlusExt}.}
  \item{adv_span}{Unit: kb, only used when rule is \code{basalPlusExt}.}
  \item{adv_twoDistance}{Unit: kb, only used when rule is \code{twoClosest}.}
  \item{adv_oneDistance}{Unit: kb, only used when rule is \code{oneClosest}.}
  \item{request_interval}{Time interval for two requests. Default is 300 seconds.}
  \item{max_tries}{Maximal times for aotumatically reconnecting GREAT web server.}
  \item{version}{Version of GREAT. The value should be "4.0.4", "3.0.0", "2.0.2". Shorten version numbers can also be used, such as using "4" or "4.0" is same as "4.0.4".}
  \item{base_url}{the url of \code{cgi-bin} path, only used when it is explicitly specified.}
  \item{help}{Whether to print help messages.}

}
\details{
Note: On Aug 19 2019 GREAT released version 4(\url{https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655442/Version+History} ) where it supports \code{hg38} genome and removes some ontologies such pathways. \code{\link{submitGreatJob}} still
takes \code{hg19} as default. \code{hg38} can be specified by the \code{species = "hg38"} argument.
To use the older versions such as 3.0.0, specify as \code{submitGreatJob(..., version = "3.0.0")}.

Note it does not use the standard GREAT API. This function directly send data to GREAT web server
by HTTP POST.

Following text is copied from GREAT web site ( \url{http://great.stanford.edu/public/html/} )

Explanation of \code{rule} and settings with names started with 'adv_' (advanced settings):

\describe{
  \item{basalPlusExt}{Mode 'Basal plus extension'. Gene regulatory domain definition:  Each gene is assigned a basal regulatory domain of a minimum distance upstream  and downstream of the TSS (regardless of other nearby genes, controlled by \code{adv_upstream} and  \code{adv_downstream} argument). The gene regulatory domain is extended in both directions  to the nearest gene's basal domain but no more than the maximum extension in one direction (controlled by \code{adv_span}).}
  \item{twoClosest}{Mode 'Two nearest genes'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the nearest  gene's TSS (controlled by \code{adv_twoDistance}) but no more than the maximum extension in one direction.}
  \item{oneClosest}{Mode 'Single nearest gene'. Gene regulatory domain definition:  Each gene is assigned a regulatory domain that extends in both directions to the midpoint  between the gene's TSS and the nearest gene's TSS (controlled by \code{adv_oneDistance}) but no more than the maximum  extension in one direction.}
}
}
\section{When_background_regions_are_set}{
Note when \code{bg} argument is set to a list of background regions, GREAT uses a completely different test!

When \code{bg} is set, \code{gr} should be exactly subset of \code{bg}. For example, let's say a background region list contains
five regions: \code{[1, 10], [15, 23], [34, 38], [40, 49], [54, 63]}, \code{gr} can only be a subset of the five regions, which
means \code{gr} can take \code{[15, 23], [40, 49]}, but it cannot take \code{[16, 20], [39, 51]}. In this setting, regions are taken
as single units and Fisher's exact test is applied for calculating the enrichment (by testing number of regions in the 2x2 contigency table).

Check \url{https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655452/File+Formats#FileFormats-Whatshouldmybackgroundregionsfilecontain?} for more explanations.

Please note from rGREAT 1.99.0, setting \code{bg} is not supported any more and this argument will be removed in the future. You can either directly use GREAT website or use other Bioconductor packages such as "LOLA" to perform
the Fisher's exact test-based analysis.

If you want to restrict the input regions to background regions (by intersections) and still to apply Binomial test there, please
consider to use local GREAT by \code{\link{great}}.}
\value{
A \code{\link{GreatJob-class}} object which can be used to get results from GREAT server. The following methods can be applied on it:

\itemize{
  \item \code{\link{getEnrichmentTables,GreatObject-method}} to retreive the result tables. 
  \item \code{\link{getRegionGeneAssociations,GreatObject-method}} to get the associations between input regions and genes.
  \item \code{\link{plotRegionGeneAssociations,GreatObject-method}} to plot the associations bewteen input regions and genes.
  \item \code{\link{shinyReport,GreatObject-method}} to view the results by a shiny application.
}
}
\seealso{
\code{\link{great}} for the local implementation of GREAT algorithm.
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
gr = randomRegions(nr = 1000, genome = "hg19")
job = submitGreatJob(gr)
job

# more parameters can be set for the job
if(FALSE) { # suppress running it when building the package
    # current GREAT version is 4.0.4
    job = submitGreatJob(gr, genome = "hg19")
    job = submitGreatJob(gr, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
    job = submitGreatJob(gr, rule = "twoClosest", adv_twoDistance = 2000)
    job = submitGreatJob(gr, rule = "oneClosest", adv_oneDistance = 2000)
}
}
