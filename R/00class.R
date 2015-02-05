## reference class for GreatJob


# == title
# Class to store and retrieve GREAT results
#
# == details
# After submitting request to GREAT server, the generated results will be 
# available on GREAT server for some time. The \code{GreatJob} class is defined 
# to store parameters that user has set and result tables what are retrieved 
# from GREAT server.
#
# == contructor
# `submitGreatJob` is used to generate a ``GreatJob`` instance.
#
GreatJob = setClass("GreatJob",
    slots = list(
        job_env = "environment",
        parameters = "list",   # parameters that are sent to GREAT
        enrichment_tables = "environment",
        association_tables = "environment",
        tempdir = "character"),
    prototype = list(
        job_env = new.env(),
        parameters = list(),
        enrichment_tables = new.env(),
        association_tables = new.env())
)

setMethod(f = "show",
          signature = "GreatJob",
          definition = function(object) {
    
    cat("Submit time:", 
        format(object@job_env$submit_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Session ID:", id(object), "\n")
    cat("Species:", param(object, "species"), "\n")
    cat("Background:", param(object, "bgChoice"), "\n")
    if(param(object, "rule") == "basalPlusExt") {
        cat("Model:", "Basal plus extension", "\n")
        cat("  Proximal:", param(object, "adv_upstream"), "kb upstream,", 
            param(object, "adv_downstream"), "kb downstream,\n  plus Distal: up to", 
            param(object, "adv_span"), "kb\n")
    } else if(param(object, "rule") == "twoClosest") {
        cat("Model:", "Two nearest genes", "\n")
        cat("  within", param(object, "adv_twoDistance"), "kb\n")
    } else if(param(object, "rule") == "oneClosest") {
        cat("Model:", "Single nearest gene", "\n")
        cat("  within", param(object, "adv_oneDistance"), "kb\n")
    }
    if(param(object, "includeCuratedRegDoms")) {
        cat("Include curated regulatory domains\n")
    }
    cat("\n")
    cat("Enrichment tables for following ontologies have been downloaded:\n")
    if(length(ls(envir = object@enrichment_tables)) == 0) {
        cat("  None\n")
    } else {
        for(nm in ls( envir = object@enrichment_tables)) {
            cat("  ", nm, "\n", sep = "")
        }
    }
    cat("\n")
})

setGeneric(name = "id",
    def = function(job) {
        standardGeneric("id")
})

setMethod(f = "id",
    signature = "GreatJob",
    definition = function(job) {
        job@job_env$id
})

setGeneric(name = "param",
    def = function(job, name) {
        standardGeneric("param")
})

setMethod(f = "param",
    signature = "GreatJob",
    definition = function(job, name) {
        job@parameters[[name]]
})
