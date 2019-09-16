
.onAttach = function(libname, pkgname) {
	msg = 
"Note: On Aug 19 2019 GREAT released version 4 where it supports `hg38` genome and removes some ontologies such pathways. `submitGreatJob()` still
takes `hg19` as default. `hg38` can be specified by the `species = 'hg38'` argument.
To use the older versions such as 3.0.0, specify as `submitGreatJob(..., version = '3.0.0')`."
    msg = paste(strwrap(msg), collapse = "\n")
    msg = paste0("\n------------------\n", msg, "\n------------------")
    packageStartupMessage(msg)
}
