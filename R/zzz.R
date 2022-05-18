
.onAttach = function(libname, pkgname) {
	msg = 
"Note: On Aug 19 2019 GREAT released version 4 where it supports `hg38` genome and removes some ontologies such pathways. `submitGreatJob()` still
takes `hg19` as default. `hg38` can be specified by the `species = 'hg38'` argument.
To use the older versions such as 3.0.0, specify as `submitGreatJob(..., version = '3.0.0')`."
    msg = paste(strwrap(msg), collapse = "\n")

    msg2 = "From rGREAT version 1.99.0, it also implements the GREAT algorithm and it allows to integrate GREAT analysis with the Bioconductor annotation ecosystem. Check the new function `great()` and the new vignette."
    
    msg2 = paste(strwrap(msg2), collapse = "\n")

    msg = paste0("\n------------------\n", msg, "\n\n", msg2, "\n------------------")

    packageStartupMessage(msg)
}
