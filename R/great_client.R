


# == title
# Class to store and retrieve GREAT results
#
# == details
# After submitting request to GREAT server, the generated results will be 
# available on GREAT server for some time. The ``GreatJob-class`` is defined 
# to store parameters that user has set and result tables what were retrieved 
# from GREAT server.
#
# == Constructor
# Users don't need to construct by hand, `submitGreatJob` is used to generate a ``GreatJob-class`` instance.
#
# == Workflow
# After submitting request to GREAT server, users can perform following steps:
#
# - `getEnrichmentTables,GreatJob-method` to get enrichment tables for selected ontologies catalogues.
# - `plotRegionGeneAssociations,GreatJob-method` to plot associations between regions and genes
# - `getRegionGeneAssociations,GreatJob-method` to get a `GenomicRanges::GRanges` object which contains associations bewteen regions and genes.
# - `shinyReport,GreatJob-method` to view the results by a shiny application.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
GreatJob = setClass("GreatJob",
    slots = list(
        job_env = "environment",
        parameters = "list",   # parameters that are sent to GREAT
        enrichment_tables = "environment",
        association_tables = "environment")
)

# == title
# Constructor method for GreatJob class
#
# == param
# -... arguments.
#
# == details
# There is no public constructor method for the `GreatJob-class`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
GreatJob = function(...) {
    new("GreatJob", ...)
}


# == title
# Perform online GREAT analysis
#
# == param
# -gr A `GenomicRanges::GRanges` object or a data frame which contains at least three columns (chr, start and end).
# -bg Not supported any more. See explanations in section "When_background_regions_are_set".
# -gr_is_zero_based Are start positions in ``gr`` zero-based?
# -species Species. "hg38", "hg19", "mm10", "mm9" are supported in GREAT version 4.x.x, "hg19", "mm10", "mm9", "danRer7" are supported in GREAT version 3.x.x and "hg19", "hg18", "mm9", "danRer7" are supported in GREAT version 2.x.x.
# -includeCuratedRegDoms  Whether to include curated regulatory domains, see https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules#AssociationRules-CuratedRegulatoryDomains .
# -rule How to associate genomic regions to genes. See 'Details' section.
# -adv_upstream Unit: kb, only used when rule is ``basalPlusExt``.
# -adv_downstream Unit: kb, only used when rule is ``basalPlusExt``.
# -adv_span Unit: kb, only used when rule is ``basalPlusExt``.
# -adv_twoDistance Unit: kb, only used when rule is ``twoClosest``.
# -adv_oneDistance Unit: kb, only used when rule is ``oneClosest``.
# -request_interval Time interval for two requests. Default is 300 seconds.
# -max_tries Maximal times for aotumatically reconnecting GREAT web server.
# -version Version of GREAT. The value should be "4.0.4", "3.0.0", "2.0.2". Shorten version numbers
#          can also be used, such as using "4" or "4.0" is same as "4.0.4".
# -base_url the url of ``cgi-bin`` path, only used when it is explicitly specified.
# -help Whether to print help messages.
#
# == details
# Note: On Aug 19 2019 GREAT released version 4(https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655442/Version+History ) where it supports ``hg38`` genome and removes some ontologies such pathways. `submitGreatJob` still
# takes ``hg19`` as default. ``hg38`` can be specified by the ``species = "hg38"`` argument.
# To use the older versions such as 3.0.0, specify as ``submitGreatJob(..., version = "3.0.0")``.
#
# Note it does not use the standard GREAT API. This function directly send data to GREAT web server
# by HTTP POST.
#
# Following text is copied from GREAT web site ( http://great.stanford.edu/public/html/ )
#
# Explanation of ``rule`` and settings with names started with 'adv_' (advanced settings):
#
# -basalPlusExt Mode 'Basal plus extension'. Gene regulatory domain definition: 
#   Each gene is assigned a basal regulatory domain of a minimum distance upstream 
#   and downstream of the TSS (regardless of other nearby genes, controlled by ``adv_upstream`` and 
#   ``adv_downstream`` argument). The gene regulatory domain is extended in both directions 
#   to the nearest gene's basal domain but no more than the maximum extension in one direction
#   (controlled by ``adv_span``).
# -twoClosest Mode 'Two nearest genes'. Gene regulatory domain definition: 
#   Each gene is assigned a regulatory domain that extends in both directions to the nearest 
#   gene's TSS (controlled by ``adv_twoDistance``) but no more than the maximum extension in one direction.
# -oneClosest Mode 'Single nearest gene'. Gene regulatory domain definition: 
#   Each gene is assigned a regulatory domain that extends in both directions to the midpoint 
#   between the gene's TSS and the nearest gene's TSS (controlled by ``adv_oneDistance``) but no more than the maximum 
#   extension in one direction.
#
# == When_background_regions_are_set
# Note when ``bg`` argument is set to a list of background regions, GREAT uses a completely different test!
#
# When ``bg`` is set, ``gr`` should be exactly subset of ``bg``. For example, let's say a background region list contains
# five regions: ``[1, 10], [15, 23], [34, 38], [40, 49], [54, 63]``, ``gr`` can only be a subset of the five regions, which
# means ``gr`` can take ``[15, 23], [40, 49]``, but it cannot take ``[16, 20], [39, 51]``. In this setting, regions are taken
# as single units and Fisher's exact test is applied for calculating the enrichment (by testing number of regions in the 2x2 contigency table).
#
# Check https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655452/File+Formats#FileFormats-Whatshouldmybackgroundregionsfilecontain? for more explanations.
#
# Please note from rGREAT 1.99.0, setting ``bg`` is not supported any more and this argument will be removed in the future. You can either directly use GREAT website or use other Bioconductor packages such as "LOLA" to perform
# the Fisher's exact test-based analysis.
#
# If you want to restrict the input regions to background regions (by intersections) and still to apply Binomial test there, please
# consider to use local GREAT by `great`.
#
# == value
# A `GreatJob-class` object which can be used to get results from GREAT server. The following methods can be applied on it:
# 
# - `getEnrichmentTables,GreatObject-method` to retreive the result tables. 
# - `getRegionGeneAssociations,GreatObject-method` to get the associations between input regions and genes.
# - `plotRegionGeneAssociations,GreatObject-method` to plot the associations bewteen input regions and genes.
# - `shinyReport,GreatObject-method` to view the results by a shiny application.
#
# == seealso
# `great` for the local implementation of GREAT algorithm.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# gr = randomRegions(nr = 1000, genome = "hg19")
# job = submitGreatJob(gr)
# job
#
# # more parameters can be set for the job
# if(FALSE) { # suppress running it when building the package
#     # current GREAT version is 4.0.4
#     job = submitGreatJob(gr, genome = "hg19")
#     job = submitGreatJob(gr, adv_upstream = 10, adv_downstream = 2, adv_span = 2000)
#     job = submitGreatJob(gr, rule = "twoClosest", adv_twoDistance = 2000)
#     job = submitGreatJob(gr, rule = "oneClosest", adv_oneDistance = 2000)
# }
#
submitGreatJob = function(gr, bg = NULL,
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
    help = TRUE
    ) {
        
    version = version[1]
    if(!version %in% names(SPECIES)) {
        stop(paste0("'version' should be in ", paste(names(SPECIES), collapse = ", "), "."))
    }
    species  = species[1]
    if(!species %in% SPECIES[[version]]) {
        stop(paste0("GREAT with version '", version, "' only supports following species:\n  ", paste(SPECIES[[version]], collapse = ", ")))
    }
    rule = match.arg(rule)[1]

    if(help) {
        cl = as.list(match.call())
        if(version == "4.0.4" && (!"species" %in% names(cl))) {
            message_wrap('Note: On Aug 19 2019 GREAT released version 4 which supports hg38 genome and removes some ontologies such
pathways. submitGreatJob() still takes hg19 as default. hg38 can be specified by argument `species = "hg38"`. To use
the older versions such as 3.0.0, specify as submitGreatJob(..., version = "3"). Set argument `help` to `FALSE`
to turn off this message.')
        }
    }

    includeCuratedRegDoms = as.numeric(includeCuratedRegDoms[1])
    
    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")
    
    if(inherits(gr, "data.frame")) {
        df = gr
        gr = GRanges(seqnames = gr[[1]],
                     ranges = IRanges(start = gr[[2]],
                                       end = gr[[3]]))
    } else if(inherits(gr, "character")) {
        df = read.table(gr, stringsAsFactors = FALSE)
        gr = GRanges(seqnames = df[[1]],
                     ranges = IRanges(start = df[[2]],
                                       end = df[[3]]))
    }
    mcols(gr) = NULL

    bgChoice = ifelse(is.null(bg), "wholeGenome", "file")

    if(!is.null(bg)) {

        warning_wrap("From rGREAT 1.99.0, `bg` will not be supported any more because GREAT requires a special format for `gr` and `bg` if both are set, and it uses a completely different method for the enrichment analysis. Please see the documentation of `submitGreatJob()` for more explanations.")
        
        if(inherits(bg, "data.frame")) {
            df_bg = bg
            bg = GRanges(seqnames = bg[[1]],
                         ranges = IRanges(start = bg[[2]],
                                           end = bg[[3]]))
            i_strand = sapply(df_bg, function(x) all(x %in% c("+", "-", "*")))
            if(any(i_strand)) {
                strand(df_bg) = df_bg[, i_strand]
            }
        } else if(inherits(bg, "character")) {
            df_bg = read.table(bg, stringsAsFactors = FALSE)
            bg = GRanges(seqnames = df_bg[[1]],
                         ranges = IRanges(start = df_bg[[2]],
                                           end = df_bg[[3]]))
            i_strand = sapply(df_bg, function(x) all(x %in% c("+", "-", "*")))
            if(any(i_strand)) {
                strand(df_bg) = df_bg[, i_strand]
            }
        }
        mcols(bg) = NULL
    }

    # check seqnames should have 'chr' prefix
    if(species %in% c("hg19", "hg18", "hg38", "mm10", "mm9")) {
        if(!all(grepl("^chr", seqnames(gr)))) {
            stop("Chromosome names (in `gr`) should have 'chr' prefix.\n")
        }
        if(!is.null(bg)) {
            if(!all(grepl("^chr", seqnames(bg)))) {
                stop("Chromosome names (in `bg`) should have 'chr' prefix.\n")
            }
        }
    } else if(species == "danRer7") {
        if(!all(grepl("^(chr|Zv9)", seqnames(gr)))) {
            stop("Chromosome names (in `gr`) should have 'chr/Zv9' prefix.\n")
        }
        if(!is.null(bg)) {
            if(!all(grepl("^(chr|Zv9)", seqnames(bg)))) {
                stop("Chromosome names (in `bg`) should have 'chr/Zv9' prefix.\n")
            }
        }
    }

    # transform GRanges to data frame
    bed = as.data.frame(gr)
    if(!gr_is_zero_based) {
        bed[, 2] = bed[, 2] - 1
    }
    bed_bg = NULL
    if(!is.null(bg)) {
        bed_bg = as.data.frame(bg)[, 1:3]
        if(!gr_is_zero_based) {
            bed_bg[, 2] = bed_bg[, 2] - 1
        }
    }

    # check request frequency
    time_interval = as.numeric(Sys.time()) - as.numeric(rGREAT_env$LAST_REQUEST_TIME)
    if(time_interval < request_interval) {
        message(qq("Don't make too frequent requests. The time break is @{request_interval}s.\nPlease wait for @{round(request_interval - time_interval)}s for the next request.\nThe time break can be set by `request_interval` argument.\n"))
        sleep(request_interval - time_interval)
    }

    #message("sending request to GREAT web server...")
    if(missing(base_url)) {
        BASE_URL = BASE_URL_LIST[version]
    } else {
        BASE_URL = base_url
    }

    # save file into temporary files
    bed = cbind(bed[, 1:3], name = paste(bed[[1]], ":", bed[[2]], "-", bed[[3]], sep = ""))
    f_bed = tempfile(fileext = ".gz")
    con = gzcon(file(f_bed, "wb"))
    write.table(bed, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
    f_bed_bg = ""
    if(!is.null(bed_bg)) {
        bed_bg = cbind(bed_bg[, 1:3], name = paste(bed_bg[[1]], ":", bed_bg[[2]], "-", bed_bg[[3]], sep = ""))
        f_bed_bg = tempfile(fileext = ".gz")
        con = gzcon(file(f_bed_bg, "wb"))
        write.table(bed_bg, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        close(con)
    }

    i_try = 0
    while(1) {
        
        if(bgChoice == "wholeGenome") {
            error = try(response <- postForm(qq("@{BASE_URL}/greatWeb.php"),
                    "species"               = species,
                    "rule"                  = rule,
                    "span"                  = adv_span,
                    "upstream"              = adv_upstream,
                    "downstream"            = adv_downstream,
                    "twoDistance"           = adv_twoDistance,
                    "oneDistance"           = adv_oneDistance,
                    "includeCuratedRegDoms" = includeCuratedRegDoms,
                    "fgChoice"              = "file",
                    "fgFile"                = fileUpload(f_bed),
                    "bgChoice"              = bgChoice,
                    "adv_upstream"          = adv_upstream,
                    "adv_downstream"        = adv_downstream,
                    "adv_span"              = adv_span,
                    "adv_twoDistance"       = adv_twoDistance,
                    "adv_oneDistance"       = adv_oneDistance
                    ))
        } else {
            error = try(response <- postForm(qq("@{BASE_URL}/greatWeb.php"),
                    "species"               = species,
                    "rule"                  = rule,
                    "span"                  = adv_span,
                    "upstream"              = adv_upstream,
                    "downstream"            = adv_downstream,
                    "twoDistance"           = adv_twoDistance,
                    "oneDistance"           = adv_oneDistance,
                    "includeCuratedRegDoms" = includeCuratedRegDoms,
                    "fgChoice"              = "file",
                    "fgFile"                = fileUpload(f_bed),
                    "bgChoice"              = "file",
                    "bgFile"                = fileUpload(f_bed_bg),
                    "adv_upstream"          = adv_upstream,
                    "adv_downstream"        = adv_downstream,
                    "adv_span"              = adv_span,
                    "adv_twoDistance"       = adv_twoDistance,
                    "adv_oneDistance"       = adv_oneDistance
                    ))
        }
        
        i_try = i_try + 1
        
        if(!inherits(error, "try-error")) {
            break
        } else {
            message(error, appendLF = FALSE)
            if(i_try > max_tries) {
                stop(qq("max try: @{max_tries} reached. Stop with error.\n"))
            } else {
                message(qq("failed with the request, try after @{ceiling(request_interval/60)} min (@{i_try}@{ifelse(i_try %% 10 == 1, 'st', ifelse(i_try %% 10 == 2, 'nd', 'th'))} try)"))
                sleep(request_interval)
            }
        }
    }

    file.remove(f_bed)
    if(!is.null(bg)) {
        file.remove(f_bed_bg)
    }

    rGREAT_env$LAST_REQUEST_TIME = Sys.time()

    # parsing error
    if(any(grepl("encountered a user error", response))) {
        if(grepl("^.*<blockquote><b>Error details:<\\/b>(.*?)<\\/blockquote>.*$", response)) {
            msg = gsub("^.*<blockquote><b>Error details:<\\/b>(.*?)<\\/blockquote>.*$", "\\1", response)
        } else {
            msg = gsub("^.*<blockquote>(.*?)<\\/blockquote>.*$", "\\1", response)
        }
        msg = gsub("<.*?>", "", msg)
        msg = gsub(" +", " ", msg)
        msg = strwrap(msg)
        msg = paste(msg, collapse = "\n")

        if(grepl("The foreground set is not a subset of the background set", msg)) {
            msg = paste0(msg, "\nVisit http://great.stanford.edu/help/display/GREAT/File+Formats#FileFormats-Whatshouldmybackgroundregionsfilecontain%3F for more explanation.\n")
        }

        stop(paste0("GREAT encountered a user error (message from GREAT web server):\n", msg))
    }

    # parsing warning
    if(any(grepl("<strong>Warning:<\\/strong>", response))) {
        msg = gsub("^.*<strong>Warning:<\\/strong>(.*?)<\\/p>.*$", "\\1", response)
        msg = gsub("<.*?>", "", msg)
        msg = gsub(" +", " ", msg)

        msg = strwrap(msg)
        msg = paste(msg, collapse = "\n")

        warning(paste0("GREAT gives a warning:\n", msg))
    }
    
    jobid = gsub("^.*var _sessionName = \"(.*?)\";.*$", "\\1",  response)
    jobid = as.vector(jobid)
      
    job = GreatJob()
    job@job_env = new.env()
    job@enrichment_tables = new.env()
    job@association_tables = new.env()

    if(version == "default") version = DEFAULT_VERSION
    
    job@parameters = list(
            "species"               = species,
            "rule"                  = rule,
            "span"                  = adv_span,
            "upstream"              = adv_upstream,
            "downstream"            = adv_downstream,
            "twoDistance"           = adv_twoDistance,
            "oneDistance"           = adv_oneDistance,
            "includeCuratedRegDoms" = includeCuratedRegDoms,
            "fgChoice"              = "data",
            "bgChoice"              = bgChoice,
            "adv_upstream"          = adv_upstream,
            "adv_downstream"        = adv_downstream,
            "adv_span"              = adv_span,
            "adv_twoDistance"       = adv_twoDistance,
            "adv_oneDistance"       = adv_oneDistance,
            "n_region"              = nrow(bed),
            "version"               = version,
            "f_bed"                 = f_bed,
            "f_bed_bg"              = f_bed_bg)
    job@parameters$n_bg = nrow(bed_bg)
    
    job@job_env$id = jobid
    job@job_env$submit_time = Sys.time()
    job@job_env$gr = bed
    job@job_env$bg = bed_bg

    td = tempdir()
    dir.create(td, showWarnings = FALSE)
    job@job_env$tempdir = td

    ### parse the response html code and extract the ontogolies
    match = gregexpr("<div class=\"gc_header gc_column\">.*?<\\/div>", response)[[1]]
    match.length = attr(match, "match.length")

    category = sapply(seq_along(match), function(i) {
            substr(response, match[i], match[i]+match.length[i]-1)
    })

    category = gsub("<.*?>", "", category)

    match = gregexpr("<div class=\"gc_info_item gc_column\">.*?<\\/div>", response)[[1]]
    match.length = attr(match, "match.length")

    term = sapply(seq_along(match), function(i) {
            substr(response, match[i], match[i]+match.length[i]-1)
    })

    term_value = gsub("^.*value=\"(\\w+?)\".*$", "\\1", term)
    term_value = gsub("&nbsp;", "", term_value)
    term_value = gsub("<.*?>", "", term_value)
    term_value = gsub("^\\s+|\\s+$", "", term_value)


    term_name = gsub("<.*?>", "", term)
    term_name = gsub("&nbsp;", "", term_name)
    term_name = gsub("^\\s+|\\s+$", "", term_name)

    ONTOLOGY_KEYS = term_value
    names(ONTOLOGY_KEYS) = term_name
    ONTOLOGY_KEYS = ONTOLOGY_KEYS[ONTOLOGY_KEYS != ""]

    CATEGORY = lapply(seq_along(category), function(i) {
            index = seq_along(term_name)
            tn = term_name[index %% length(category) == i %% length(category)]
            tn[tn != ""]
    })
    names(CATEGORY) = category

    job@job_env$CATEGORY = CATEGORY
    job@job_env$ONTOLOGY_KEYS = ONTOLOGY_KEYS

    return(job)
}

## reference class for GreatJob


setMethod(f = "show",
          signature = "GreatJob",
          definition = function(object) {
    
    cat("Submit time:", 
        format(object@job_env$submit_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("  Note the results may only be avaiable on GREAT server for 24 hours.\n")
    cat("Version:", param(object, "version"), "\n")
    #cat("Session ID:", id(object), "\n")
    cat("Species:", param(object, "species"), "\n")
    cat("Inputs:", param(object, "n_region"), "regions\n")
    # if(param(object, "bgChoice") == "wholeGenome") {
    #     cat("Background:", param(object, "bgChoice"), "\n")
    # } else {
    #     cat("Background: user-defined,", param(object, "n_bg"), "regions\n")
    # }
    if(param(object, "rule") == "basalPlusExt") {
        cat("Mode:", "Basal plus extension", "\n")
        cat("  Proximal:", param(object, "adv_upstream"), "kb upstream,", 
            param(object, "adv_downstream"), "kb downstream,\n  plus Distal: up to", 
            param(object, "adv_span"), "kb\n")
    } else if(param(object, "rule") == "twoClosest") {
        cat("Mode:", "Two nearest genes", "\n")
        cat("  within", param(object, "adv_twoDistance"), "kb\n")
    } else if(param(object, "rule") == "oneClosest") {
        cat("Mode:", "Single nearest gene", "\n")
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

# == title
# Get enrichment tables from GREAT web server
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -ontology Ontology names. Valid values are in `availableOntologies`. ``ontology`` is prior to 
#           ``category`` argument.
# -category Pre-defined ontology categories. One category can contain more than one ontologies. Valid values are in 
#            `availableCategories`
# -request_interval Time interval for two requests. Default is 300 seconds.
# -max_tries Maximal times for automatically reconnecting GREAT web server.
# -download_by Internally used.
# -verbose Whether to print messages.
#
# == value
# The structure of the data frames are same as the tables available on GREAT website.
#
# == see also
# `availableOntologies`, `availableCategories`
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
# tbl = getEnrichmentTables(job)
# names(tbl)
# head(tbl[[1]])
# job
#
# tbl = getEnrichmentTables(job, ontology = "GO Molecular Function")
# tbl = getEnrichmentTables(job, category = "GO")
#
setMethod(f = "getEnrichmentTables",
    signature = "GreatJob",
    definition = function(object, ontology = NULL, category = "GO",
    request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
    verbose = TRUE) {
        
    job = object

    jobid = id(job)
    species = param(job, "species")
    version = param(job, "version")
    BASE_URL = BASE_URL_LIST[version]

    if(!file.exists(job@job_env$tempdir)) {
        td = tempdir()
        dir.create(td, showWarnings = FALSE)
        job@job_env$tempdir = td
    }

    if(is.null(ontology)) {
        if(is.null(category)) {
            stop("`ontology` and `category` can not be both NULL.\n")
        }
        if(length(setdiff(category, availableCategories(job)))) {
            stop("Only categories in `availableCategories()` are allowed.\n")
        }
        ontology = availableOntologies(job, category = category)
    }
    ontology = unique(ontology)
    
    if(length(setdiff(ontology, availableOntologies(job)))) {
        stop("Only ontologies in `availableOntologies()` are allowed.\n")
    }

    ONTOLOGY_KEYS = job@job_env$ONTOLOGY_KEYS
    
    download_by = match.arg(download_by)[1]

    if(download_by == "json" && verbose) {
        message_wrap("The default enrichment tables contain no associated genes for the input regions.",
            "You can set `download_by = 'tsv'` to download the complete table,",
            "but note only the top 500 regions can be retreived. See the following link:\n\n",
            "https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655401/Export#Export-GlobalExport\n")
    }

    res = lapply(ontology, function(onto) {
        if(!is.null(job@enrichment_tables[[onto]]) && download_by == "json") {
            res = job@enrichment_tables[[onto]]
            return(res)
        } else {
            if(download_by == "json") {
                res = GREAT.read.json(job, qq("@{BASE_URL}/readJsFromFile.php?path=/scratch/great/tmp/results/@{jobid}.d/@{ONTOLOGY_KEYS[onto]}.js"), onto, 
                    request_interval = request_interval, max_tries = max_tries)
                job@enrichment_tables[[onto]] = res
            } else {
                res = download_enrichment_table(job, ONTOLOGY_KEYS[onto], request_interval = request_interval, max_tries = max_tries)
            }
            return(res)
        }
    })
    names(res) = ontology
    return(res)
})



# == title
# Get a single enrichment table from GREAT web server
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -ontology A single ontology names. Valid values are in `availableOntologies`. 
# -... All pass to `getEnrichmentTables,GreatJob-method`.
#
# == value
# A data frame of the enrichment results for a single ontology.
# 
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
# tb = getEnrichmentTable(job, ontology = "GO Molecular Function")
# head(tb)
setMethod(f = "getEnrichmentTable",
    signature = "GreatJob",
    definition = function(object, ontology, ...) {

    if(length(ontology) != 1) {
        stop_wrap("Length of `ontology` must be 1.")
    }
    getEnrichmentTables(object, ontology, ...)[[1]]
})

message_wrap = function (...)  {
    x = paste(..., collapse = " ")
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}

download_enrichment_table = function(job, onto, request_interval = 10, max_tries = 100) {
    jobid = id(job)
    species = param(job, "species")
    bgChoice = param(job, "bgChoice")
    version = param(job, "version")
    BASE_URL = BASE_URL_LIST[version]

    i_try = 0
    while(1) {
        error = try(
            response <- postForm(qq("@{BASE_URL}/downloadAllTSV.php"),
                outputDir = qq("/scratch/great/tmp/results/@{jobid}.d/"),
                outputPath = qq("/tmp/results/@{jobid}.d/"),
                ontoName = onto,
                species = species,
                ontoList = qq("@{onto}@1-Inf"),
                binom = ifelse(bgChoice == "wholeGenome", "true", "false")
            )
        )
            
        i_try = i_try + 1
            
        if(!inherits(error, "try-error")) {
            break
        } else {
                
            if(i_try > max_tries) {
                stop(qq("max try: @{max_tries} reached. Stop with an error.\n"))
            } else {
                message("failed to download, try after 30s")
                sleep(request_interval)
            }
        }
    }

    response = strsplit(response, "\n")[[1]]
    response = response[-(1:3)]
    response[[1]] = gsub("^#", "", response[1])
    response = response[!grepl("^#", response)]
    response = paste(response, collapse = "\n")
    error = try(tb <- read.table(textConnection(response), sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE))
    if(inherits(error, "try-error")) {
        stop("downloading enrichment table failed.")
    }

    return(tb)
}

# == title
# All available ontology names of the GREAT job
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -category one or multiple categories. All available categories can be got by `availableCategories`.
#
# == details
# The values of the supported ontologies sometime change. You should run the function to get the real-time
# values. The meaning of ontology returned is quite self-explained by the name.
#
# == value
# The returned values is a vector of ontologies.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
# availableOntologies(job)
# availableOntologies(job, category = "Pathway Data")
#
setMethod(f = "availableOntologies",
    signature = "GreatJob",
    definition = function(object, category = NULL) {
        
    job = object

    species = param(job, "species")
    CATEGORY = job@job_env$CATEGORY
    
    if(is.null(category)) {
        onto = unlist(CATEGORY)
        names(onto) = NULL
        return(onto)
    } else {
        if(length(setdiff(category, availableCategories(job))) > 0) {
            stop("Value of `category` is invalid. Please use availableCategories(job) to find supported categories.\n")
        }
        onto = unlist(CATEGORY[category])
        names(onto) = NULL
        return(onto)
    }
})


# == title
# Available ontology categories of the GREAT job
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
#
# == details
# The values of the supported categories sometime change. You should run the function to get the real-time
# values. The meaning of categories returned is quite self-explained by the name.
#
# == value
# The returned value is a vector of categories.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
# availableCategories(job)
#
setMethod(f = "availableCategories",
    signature = "GreatJob",
    definition = function(object) {

    species = param(object, "species")
    CATEGORY = object@job_env$CATEGORY
    
    names(CATEGORY)
})

# download from `url` and save to `file`
download = function(url, file, request_interval = 10, max_tries = 100) {

    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")


    #message(qq("try to download from @{url}"))
    i_try = 0
    while(1) {
        error = try(download.file(url, destfile = file, quiet = TRUE))
            
        i_try = i_try + 1
            
        if(!inherits(error, "try-error")) {
            break
        } else {
            if(file.exists(file)) file.remove(file)
                
            if(i_try > max_tries) {
                stop(qq("max try: @{max_tries} reached. Stop with an error.\n"))
            } else {
                message("failed to download, try after 30s")
                sleep(request_interval)
            }
        }
    }
}

GREAT.read.json = function(job, url, onto, request_interval = 10, max_tries = 100) {
    jobid = id(job)

    if(!file.exists(job@job_env$tempdir)) {
        td = tempdir()
        dir.create(td, showWarnings = FALSE)
        job@job_env$tempdir = td
    }

    TEMP_DIR = job@job_env$tempdir

    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")

    
    f1 = qq("@{TEMP_DIR}/@{jobid}_@{onto}.js")
    f1 = gsub(" ", "_", f1)

    download(url, file = f1, request_interval = request_interval, max_tries = max_tries)

    # if something wrong, the returned js has length of 0
    if(file.info(f1)$size == 0) {
        stop("Retrieved 0 byte from remote server. Probably your job on GREAT server expires\nand you need to submit it again.\n")
    }
    
    # just in case downloading json file is interrupted
    error = try(json <- fromJSON(file = f1))
    if(inherits(error, "try-error")) {
        stop("Downloading seems interrupted. Please re-run your command.\n")
    }
    
    file.remove(f1)
    
    if(length(json) == 0) {
        return(NULL)
    }
    lt = vector("list", length(json[[1]]))
    
    for(i in seq_along(json)) {
        for(j in seq_along(json[[i]])) {
            lt[[j]] = c(lt[[j]], json[[i]][[j]])
        }
    }
    res = as.data.frame(lt, stringsAsFactors = FALSE)

    if (param(job, "bgChoice") == "file") {
      colnames(res) = c("ID", "name", "Hyper_Total_Regions", "Hyper_Expected", "Hyper_Foreground_Region_Hits",
                        "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Term_Region_Coverage",
                        "Hyper_Foreground_Gene_Hits", "Hyper_Background_Gene_Hits", "Total_Genes_Annotated", "Hyper_Raw_PValue")
      res$Hyper_Adjp_BH = p.adjust(res$Hyper_Raw_PValue, method = "BH")
    } else {
      colnames(res) = c("ID", "name", "Binom_Genome_Fraction", "Binom_Expected", "Binom_Observed_Region_Hits", "Binom_Fold_Enrichment",
                        "Binom_Region_Set_Coverage", "Binom_Raw_PValue", "Hyper_Total_Genes", "Hyper_Expected",
                        "Hyper_Observed_Gene_Hits", "Hyper_Fold_Enrichment", "Hyper_Gene_Set_Coverage",
                        "Hyper_Term_Gene_Coverage", "Hyper_Raw_PValue")
      res$Binom_Adjp_BH = p.adjust(res$Binom_Raw_PValue, method = "BH")
      res$Hyper_Adjp_BH = p.adjust(res$Hyper_Raw_PValue, method = "BH")
      res = res[, c("ID", "name", "Binom_Genome_Fraction", "Binom_Expected", "Binom_Observed_Region_Hits", "Binom_Fold_Enrichment",
                        "Binom_Region_Set_Coverage", "Binom_Raw_PValue", "Binom_Adjp_BH", "Hyper_Total_Genes", "Hyper_Expected",
                        "Hyper_Observed_Gene_Hits", "Hyper_Fold_Enrichment", "Hyper_Gene_Set_Coverage",
                        "Hyper_Term_Gene_Coverage", "Hyper_Raw_PValue", "Hyper_Adjp_BH")]
    }
    return(res)
}

# transform the original file into the standard data frame
parseRegionGeneAssociationFile = function(f1) {

    error = try(data <- read.table(f1, sep = "\t", stringsAsFactors = FALSE))
    if(inherits(error, "try-error")) {
        stop("Downloading seems interrupted. Please re-run your command.\n")
    }
    
    region = NULL
    gene = NULL
    distance = NULL
    for(i in seq_len(nrow(data))) {
        if(data[i, 2] == "NONE") {
            region = c(region, data[i, 1])
            gene = c(gene, NA)
            distance = c(distance, NA)
        } else {
            r = gregexpr("([a-zA-Z0-9\\-_.]+) \\(([+-]?\\d+)\\)", data[i, 2], perl = TRUE)[[1]]
            capture.start = attr(r, "capture.start")
            capture.length = attr(r, "capture.length")
            k = nrow(capture.start)
            region = c(region, rep(data[i, 1], k))
            gene = c(gene, substr(rep(data[i, 2], k), capture.start[, 1], capture.start[, 1] + capture.length[, 1] - 1))
            distance = c(distance, substr(rep(data[i, 2], k), capture.start[, 2], capture.start[, 2] + capture.length[, 2] - 1))
        }
    }

    da = strsplit(region, "[-:]")
    chr = sapply(da, function(x) x[1])
    start = as.integer(sapply(da, function(x) x[2]))
    end = as.integer(sapply(da, function(x) x[3]))
    df = data.frame(chr = chr,
                    start = start,
                    end = end,
                    gene = gene,
                    distTSS = as.numeric(distance),
                    stringsAsFactors = FALSE)
    rownames(df) = NULL
    return(df)    
}


# == title
# Plot region-gene associations
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -ontology A single ontology names. Valid values are in `availableOntologies`. 
# -term_id Term id in the selected ontology
# -which_plot Which plots to draw? The value should be in ``1, 2, 3``. See "Details" section for explanation.
# -request_interval Time interval for two requests. Default is 300 seconds.
# -max_tries Maximal times for automatically reconnecting GREAT web server.
# -verbose Whether to show messages.
#
# == details
# There are following figures:  
#
# - Association between regions and genes (``which_plot = 1``).
# - Distribution of distance to TSS (``which_plot = 2``).
# - Distribution of absolute distance to TSS (``which_plot = 3``).
#
# If ``ontology`` and ``term_id`` are set, only regions and genes corresponding to 
# selected ontology term will be used. Valid value for ``ontology`` is in 
# `availableOntologies` and valid value for ``term_id`` is from 'id' column 
# in the table which is returned by `getEnrichmentTables`.  
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
#
# plotRegionGeneAssociations(job)
# plotRegionGeneAssociations(job, which_plot = 1)
# plotRegionGeneAssociations(job, ontology = "GO Molecular Function",
#     term_id = "GO:0004984")
#
setMethod(f = "plotRegionGeneAssociations",
    signature = "GreatJob",
    definition = function(object, ontology = NULL, term_id = NULL, which_plot = 1:3, 
    request_interval = 10, max_tries = 100, verbose = great_opt$verbose) {

    job = object
    termID = term_id

    gr_all = getRegionGeneAssociations(object, request_interval = request_interval,
        max_tries = max_tries, verbose = verbose)

    if(!is.null(term_id)) {
        gr_term = getRegionGeneAssociations(object, ontology = ontology, term_id = term_id, request_interval = request_interval,
            max_tries = max_tries, verbose = verbose)
    } else {
        gr_term = NULL
    }

    plot_great(gr_all, gr_term, which_plot = which_plot, gr_full_len = nrow(object@job_env$gr), term_id = term_id)
})


plot_great = function(gr_all, gr_term = NULL, which_plot = 1:3, gr_full_len, term_id = NULL) {

    op = par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, length(intersect(which_plot, 1:3))), mar = c(6, 4, 4, 1))

    using_term = !is.null(gr_term)

    df_all = data.frame(distTSS = unlist(gr_all$dist_to_TSS))
    if(using_term) {
        df_term = data.frame(distTSS = unlist(gr_term$dist_to_TSS))
    }
   
    # make plots
    if(1 %in% which_plot) {
        if(using_term) {
            tb = table(table(unlist(gr_term$annotated_genes)))
            vt = numeric(10)
            vt[as.numeric(names(tb))] = tb
            vt[is.na(vt)] = 0
            v = c(vt[1:9], sum(vt[10:length(vt)]))
            names(v) = c(as.character(1:9), ">= 10")
            v[is.na(v)] = 0
            p = v/sum(v)
            pos = barplot(p, col = "black", xlab = "Number of associated regions per gene", ylab = "This term's genes", ylim = c(0, max(p)*1.5), main = qq("Number of associated regions per gene\nTerm: @{term_id}"))
            text(pos[, 1], p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        } else {
            tb = table(table(unlist(gr_all$annotated_genes)))
            v = c(gr_full_len - length(gr_all), tb["1"], tb["2"], sum(tb[as.numeric(names(tb)) > 2]))
            names(v) = c("0", "1", "2", "> 3")
            v[is.na(v)] = 0
            p = v/sum(v)
            pos = barplot(p, col = c("red", "grey", "grey", "grey"), xlab = "Number of associated genes per region", ylab = "Genomic regions", ylim = c(0, max(p)*1.5), main = "Number of associated genes per region")
            text(pos[, 1], p + 0.01, v, adj = c(0.5, 0), col = c("red", "black", "black", "black"), cex = 0.8)
            legend("topright", pch = 15, col = c("grey", "red"), legend = c("Genomic regions associated with one or more genes", "Genomic regions not associated with any genes"), cex = 0.8)
        }
    }
    if(2 %in% which_plot) {
        v = cbind(
            c("<-500"       = sum(df_all$distTSS < -500000),
              "[-500, -50)" = sum(df_all$distTSS >= -500000 & df_all$distTSS < -50000),
              "[-50, -5)"   = sum(df_all$distTSS >= -50000  & df_all$distTSS < -5000),
              "[-5, 0]"     = sum(df_all$distTSS >= -5000   & df_all$distTSS < 0),
              "0"           = sum(df_all$distTSS == 0),
              "(0, 5]"      = sum(df_all$distTSS > 0       & df_all$distTSS <= 5000),
              "(5, 50]"     = sum(df_all$distTSS > 5000    & df_all$distTSS <= 50000),
              "(50, 500]"   = sum(df_all$distTSS > 50000   & df_all$distTSS <= 500000),
              "> 500"       = sum(df_all$distTSS > 500000)))
        
        if(using_term) {
            v = cbind( 
            c("<-500"       = sum(df_term$distTSS < -500000),
              "[-500, -50)" = sum(df_term$distTSS >= -500000 & df_term$distTSS < -50000),
              "[-50, -5)"   = sum(df_term$distTSS >= -50000  & df_term$distTSS < -5000),
              "[-5, 0]"     = sum(df_term$distTSS >= -5000   & df_term$distTSS < 0),
              "0"           = sum(df_term$distTSS == 0),
              "(0, 5]"      = sum(df_term$distTSS > 0       & df_term$distTSS <= 5000),
              "(5, 50]"     = sum(df_term$distTSS > 5000    & df_term$distTSS <= 50000),
              "(50, 500]"   = sum(df_term$distTSS > 50000   & df_term$distTSS <= 500000),
              "> 500"       = sum(df_term$distTSS > 500000)), v)
        }
        p = v
        for(ic in seq_len(ncol(v))) {
            p[, ic] = v[, ic]/sum(v[, ic])
        }
        p[is.na(p)] = 0
        rownames(p) = NULL
        if(using_term) {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "", ylab = "Region-gene associations (fraction)", ylim = c(0, max(p)*1.5), main = qq("Binned by orientation and distance to TSS\nTerm: @{term_id}"), axes = FALSE)
        } else {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "", ylab = "Region-gene associations (fraction)", ylim = c(0, max(p)*1.5), main = qq("Binned by orientation and distance to TSS"), axes = FALSE)
        }
        mtext("Distance to TSS (kb)", line = 4, side = 1, cex = 0.7)
        text(t(pos), p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        axis(side = 2)
        op = par("xpd")
        par(xpd = NA)
        text(colMeans(pos), -0.02, rownames(v), srt = 45, adj = c(1, 0.5))
        par(xpd = op)
        x1 = (mean(pos[, 5]) + mean(pos[, 5]))/2
        x2 = (mean(pos[, 6]) + mean(pos[, 6]))/2
        y = max(p)*1.2
        lines( c(x1, x1), c(0, y))
        arrows(x1, y, x2, y, angle = 15, length = 0.1, code = 2)
        text(mean(pos[, 6]), y+0.01, "TSS", adj = c(0.5, 0), cex = 0.8)
        if(using_term) {
            legend("topright", pch = 15, col = c("green", "blue"), legend = c("This term's region-gene associations", "Set-wide region-gene associations"), cex = 0.8)
        }
    }
    if(3 %in% which_plot) {
        v = cbind(
            c("0"           = sum(abs(df_all$distTSS) == 0),
              "(0, 5]"      = sum(abs(df_all$distTSS) > 0       & abs(df_all$distTSS) <= 5000),
              "(5, 50]"     = sum(abs(df_all$distTSS) > 5000    & abs(df_all$distTSS) <= 50000),
              "(50, 500]"   = sum(abs(df_all$distTSS) > 50000   & abs(df_all$distTSS) <= 500000),
              "> 500"       = sum(abs(df_all$distTSS) > 500000)))
        if(using_term) {
            v = cbind( 
            c("0"           = sum(abs(df_term$distTSS) == 0),
              "(0, 5]"      = sum(abs(df_term$distTSS) > 0       & abs(df_term$distTSS) <= 5000),
              "(5, 50]"     = sum(abs(df_term$distTSS) > 5000    & abs(df_term$distTSS) <= 50000),
              "(50, 500]"   = sum(abs(df_term$distTSS) > 50000   & abs(df_term$distTSS) <= 500000),
              "> 500"       = sum(abs(df_term$distTSS) > 500000)), v)
        }
        p = v
        for(ic in seq_len(ncol(v))) {
            p[, ic] = v[, ic]/sum(v[, ic])
        }
        p[is.na(p)] = 0
        rownames(p) = NULL
        if(using_term) {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "", ylab = "Region-gene associations (fraction)", ylim = c(0, max(p)*1.5), main = qq("Binned by absolute distance to TSS\nTerm: @{term_id}"), axes = FALSE)
        } else {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "", ylab = "Region-gene associations (fraction)", ylim = c(0, max(p)*1.5), main = qq("Binned by absolute distance to TSS"), axes = FALSE)
        }
        mtext("Absolute distance to TSS (kb)", line = 4, side = 1, cex = 0.7)
        text(t(pos), p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        axis(side = 2)
        op = par("xpd")
        par(xpd = NA)
        text(colMeans(pos), -0.02, rownames(v), srt = 45, adj = c(1, 0.5))
        par(xpd = op)
        x1 = 0.9
        y = max(p)*1.2
        arrows(x1, 0, x1, y, angle = 15, length = 0.1, code = 2)
        text(x1, y+0.01, "TSS", adj = c(0, 0), cex = 0.8)
        
        if(using_term) {
            legend("topright", pch = 15, col = c("green", "blue"), legend = c("This term's region-gene associations", "Set-wide region-gene associations"), cex = 0.8)
        }
    }
}

# == title
# Plot region-gene associations
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -... All passed to `plotRegionGeneAssociations,GreatJob-method`.
#
# == details
# This function will be removed in the future, please use `plotRegionGeneAssociations,GreatJob-method` instead.
#
setMethod(f = "plotRegionGeneAssociationGraphs",
    signature = "GreatJob",
    definition = function(object, ...) {

    plotRegionGeneAssociations(object, ...)
})

# == title
# Get region-gene associations
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -ontology ontology name
# -term_id Term id in the selected ontology.
# -request_interval Time interval for two requests. Default is 300 seconds.
# -max_tries Maximal times for automatically reconnecting GREAT web server.
# -verbose Whether to show messages.
#
# == value
# A `GenomicRanges::GRanges` object. Please the two meta columns are in formats of ``CharacterList``
# and ``IntegerList`` because a region may associate to multiple genes.
#
# Please note, the distance is from the middle points of input regions to TSS.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# # note the `job` was generated from GREAT 3.0.0
# job = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
#
# gr = getRegionGeneAssociations(job)
# gr
setMethod(f = "getRegionGeneAssociations",
    signature = "GreatJob",
    definition = function(object, ontology = NULL, term_id = NULL, 
    request_interval = 10, max_tries = 100, verbose = great_opt$verbose) {

    job = object
    termID = term_id

    if(!file.exists(job@job_env$tempdir)) {
        td = tempdir()
        dir.create(td, showWarnings = FALSE)
        job@job_env$tempdir = td
    }
    
    jobid = id(job)
    species = param(job, "species")
    TEMP_DIR = job@job_env$tempdir
    version = param(job, "version")
    BASE_URL = BASE_URL_LIST[version]

    opqq = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(opqq))
    qq.options(code.pattern = "@\\{CODE\\}")

    # make plot and return value for a single term
    using_term = FALSE
    
    if(!is.null(ontology) && !is.null(termID)) {
        using_term = TRUE
        
        ontology = ontology[1]
        termID = termID[1]

        if(! ontology %in% availableOntologies(job)) {
            stop("Value of `ontology` should be in `availableOntologies()`\n")
        }
    }

    if(sum(c(is.null(ontology), is.null(termID))) == 1) {
        stop("You should set both of `ontology` and `term_id` or neither of them.\n")
    }
    
    if(using_term) {
        # check whether termID is in ontology if ontology table is already downloaded
        if(!is.null(job@enrichment_tables[[ontology]])) {
            if(! termID %in% job@enrichment_tables[[ontology]]$ID) {
                stop(qq("Cannot find '@{termID}' in enrichment table of '@{ontology}'"))
            }
        }
    }

    ONTOLOGY_KEYS = job@job_env$ONTOLOGY_KEYS

    if(using_term) {

        # prepare file names for local table
        f_term = qq("@{jobid}-@{ontology}-@{termID}-@{species}-region.txt")
        f_term = gsub('[\\/:*?"<>|]', "_", f_term)
        f_term = qq("@{TEMP_DIR}/@{f_term}")
        
        if(!is.null(job@association_tables[[qq("@{ontology}-@{termID}")]])) {
            df_term = job@association_tables[[qq("@{ontology}-@{termID}")]]
        } else {
            if(verbose) qqcat("The webpage for '@{ontology}:@{termID}' is available at:\n  @{BASE_URL}/showTermDetails.php?termId=@{termID}&ontoName=@{ONTOLOGY_KEYS[ontology]}&ontoUiName=@{ontology}&sessionName=@{jobid}&species=@{species}&foreName=@{basename(param(job, 'f_bed'))}&backName=@{basename(param(job, 'f_bed_bg'))}&table=region\n")

            if (param(job, "bgChoice") != "data") {
              url = qq("@{BASE_URL}/downloadAssociations.php?termId=@{termID}&ontoName=@{ONTOLOGY_KEYS[ontology]}&sessionName=@{jobid}&species=@{species}&foreName=@{basename(param(job, 'f_bed'))}&backName=@{basename(param(job, 'f_bed_bg'))}&table=region")
            } else {
              url = qq("@{BASE_URL}/downloadAssociations.php?termId=@{termID}&ontoName=@{ONTOLOGY_KEYS[ontology]}&sessionName=@{jobid}&species=@{species}&foreName=@{basename(param(job, 'f_bed'))}&backName=@{basename(param(job, 'f_bed_bg'))}&table=region")
            }
            download(url, file = f_term, request_interval = request_interval, max_tries = max_tries)
            check_asso_file(f_term)
            df_term = parseRegionGeneAssociationFile(f_term)
            df_term$gene[df_term$gene == ""] = NA
            job@association_tables[[qq("@{ontology}-@{termID}")]] = df_term
            file.remove(f_term)
        }
    }
    
    f_all = qq("@{TEMP_DIR}/@{jobid}-@{species}-all-region.txt")
        
    # download
    if(!is.null(job@association_tables[["all"]])) {
        df_all = job@association_tables[["all"]]
    } else {
        url = qq("@{BASE_URL}/downloadAssociations.php?sessionName=@{jobid}&species=@{species}&foreName=@{basename(param(job, 'f_bed'))}&backName=@{basename(param(job, 'f_bed_bg'))}&table=region")
        download(url, file = f_all, request_interval = request_interval, max_tries = max_tries)
        check_asso_file(f_all)
        df_all = parseRegionGeneAssociationFile(f_all)
        df_all$gene[df_all$gene == ""] = NA
        job@association_tables[["all"]] = df_all
        file.remove(f_all)
    }

    # some values of gene is NA
    if(using_term) {
        df_term_NA = df_term[is.na(df_term$gene), , drop = FALSE]
        df_term = df_term[!is.na(df_term$gene), , drop = FALSE]
    }
    df_all_NA = df_all[is.na(df_all$gene), , drop = FALSE]
    df_all = df_all[!is.na(df_all$gene), , drop = FALSE]

    
    if(using_term) {
        df = df_term
    } else {
        df = df_all
    }

    fa = paste0(df[, 1], df[, 2], df[, 3], sep = "-")
    gr = GRanges(seqnames = tapply(df[, 1], fa, function(x) x[1]),
                ranges = IRanges(tapply(df[, 2], fa, function(x) x[1]), tapply(df[, 3], fa, function(x) x[1])))
    gr$annotated_genes = vector("list", length(gr))
    gr$annotated_genes = split(df[, 4], fa)
    gr$dist_to_TSS = vector("list", length(gr))
    gr$dist_to_TSS = split(df[, 5], fa)

    gr$annotated_genes = CharacterList(gr$annotated_genes)
    gr$dist_to_TSS = IntegerList(gr$dist_to_TSS)

    seqlevels(gr) = sort_chr(unique(df[[1]]))
    gr = sort(gr)

    return(gr)
})


sleep = function(time) {
    pb = txtProgressBar(style = 3)
    for(i in seq_len(time)/time) {Sys.sleep(1); setTxtProgressBar(pb, i)}
    close(pb)
}

sort_chr = function(x) {
    y = gsub("^chr(\\d)$", "chr0\\1", x)
    y = gsub("^chr(\\d)_", "chr0\\1_", y)
    x[order(y)]
}


check_asso_file = function(file) {

    # for downloading association tables
    if(file.info(file)$size < 1000) {

        # if session expires, there is a error page
        file_text = readLines(file)
        if(any(grepl("GREAT is unable to locate a session", file_text))) {
            stop("Your job on GREAT server expires and you need to submit it again.\n")
        }
        
        # wrong termID returns a table only contain a newline
        if(all(grepl("^\\s*$", file_text))) {
            stop("Empty data, probably your 'termID' is invalid. Please check the URL above.\n")
        }
    }
}


# # == title
# # Add region-gene association to the internal enrichment tables
# #
# # == param
# # -job a `GreatJob-class` instance
# # -ontology ontology name
# # -verbose whether show message
# #
# # == details
# # After successfully run this function, users need to rerun `getEnrichmentTables` to 
# # get the enrichment tables that contains the new gene association column.
# setMethod(f = "addRegionGeneAssociation",
#     signature = "GreatJob",
#     definition = function(job, ontology = names(job@enrichment_tables), verbose = TRUE) {

#     if(length(ontology) == 0) {
#         stop("You should run `getEnrichmentTables()` first.")
#     }
#     ontology = intersect(ontology, names(job@enrichment_tables))
#     if(length(ontology) == 0) {
#         stop("Cannot find some of the ontologies.")
#     }

#     for(onto in ontology) {
#         enrichment_tb = job@enrichment_tables[[onto]]
#         all_termID = enrichment_tb$ID
#         n_term = length(all_termID)

#         genes = character(n_term)
#         for(i in seq_along(all_termID)) {
#             tid = all_termID[i]
#             if(verbose) {
#                 qqcat("Downloading associated genes for @{onto}, @{tid}, @{i}/@{n_term}...\n")
#             }
#             gr = plotRegionGeneAssociationGraphs(job, ontology = onto, termID = tid, request_interval = 5, verbose = FALSE, plot = FALSE)
#             genes[i] = paste(unique(sort(gr$gene)), collapse = ",")
#         }
#         enrichment_tb$Asso_Genes = genes
#         job@enrichment_tables[[onto]] = enrichment_tb
#     }
#     if(verbose) {
#         cat("Done. Please rerun `getEnrichmentTables()` to get the tables.\n")
#     }
#     invisible(NULL)
# })

