## reference class for GREAT_Job

GREAT_Job = setRefClass("GREAT_Job",
    fields = list(
        id = "character",
        parameters = "list",
        enrichment_tables = "list",
        association_tables = "list")
)

GREAT_Job$methods(show = function() {
    
    cat("Submit time:", 
        format(.self$parameters$submit_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Session ID:", .self$id, "\n")
    cat("Species:", .self$parameters$species, "\n")
    cat("Background:", .self$parameters$bgChoice, "\n")
    if(.self$parameters$rule == "basalPlusExt") {
        cat("Model:", "Basal plus extension", "\n")
        cat("  Proximal:", .self$parameters$adv_upstream, "kb upstream,", 
            .self$parameters$adv_downstream, "kb downstream,\n  plus Distal: up to", 
            .self$parameters$adv_span, "kb\n")
    } else if(.self$parameters$rule == "twoClosest") {
        cat("Model:", "Two nearest genes", "\n")
        cat("  within", .self$parameters$adv_twoDistance, "kb\n")
    } else if(.self$parameters$rule == "oneClosest") {
        cat("Model:", "Single nearest gene", "\n")
        cat("  within", .self$parameters$adv_oneDistance, "kb\n")
    }
    if(.self$parameters$includeCuratedRegDoms) {
        cat("Include curated regulatory domains\n")
    }
    cat("\n")
    cat("Enrichment tables for following ontologies have been downloaded:\n")
    if(length(.self$enrichment_tables) == 0) {
        cat("  None\n")
    } else {
        for(nm in names(.self$enrichment_tables)) {
            cat("  ", nm, "\n", sep = "")
        }
    }
    cat("\n")
})

GREAT_Job$methods(get_id = function() {
    .self$id
})

GREAT_Job$methods(get_param = function(name) {
    .self$parameters[[name]]
})