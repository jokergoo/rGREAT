
.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

    msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/rGREAT/
Github page: https://github.com/jokergoo/rGREAT

If you use it in published research, please cite:
Gu, Z. rGREAT: an R/Bioconductor package for functional enrichment 
  on genomic regions. Bioinformatics 2023.

This message can be suppressed by:
  suppressPackageStartupMessages(library(rGREAT))

Note: On Aug 19 2019 GREAT released version 4 where it supports `hg38`
genome and removes some ontologies such pathways. `submitGreatJob()`
still takes `hg19` as default. `hg38` can be specified by the `species
= 'hg38'` argument. To use the older versions such as 3.0.0, specify as
`submitGreatJob(..., version = '3.0.0')`.

From rGREAT version 1.99.0, it also implements the GREAT algorithm and
it allows to integrate GREAT analysis with the Bioconductor annotation
ecosystem. By default it supports more than 500 organisms. Check the
new function `great()` and the new vignette.
========================================
")  

    packageStartupMessage(msg)
}
