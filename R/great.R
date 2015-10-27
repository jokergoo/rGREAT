


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
# - call `getEnrichmentTables` to get enrichment tables for selected ontologies catalogues.
# - call `plotRegionGeneAssociationGraphs` to get associations between regions and genes
#   as well as making plots.  
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
# Send requests to GREAT web server
#
# == param
# -gr A `GenomicRanges::GRanges` object or a data frame which contains at least three columns (chr, start and end). Regions for test.
# -bg A `GenomicRanges::GRanges` object or a data frame. Background regions if needed.
# -species Species. Only four species ("hg19", "hg18", "mm9", "danRer7") are supported.
# -includeCuratedRegDoms  Whether to include curated regulatory domains.
# -bgChoice  How to define background. If it is set as ``data``, ``bg`` should be set as well.
# -rule How to associate genomic regions to genes. See 'details' section.
# -adv_upstream Unit: kb, only used when rule is ``basalPlusExt``
# -adv_downstream Unit: kb, only used when rule is ``basalPlusExt``
# -adv_span Unit: kb, only used when rule is ``basalPlusExt``
# -adv_twoDistance Unit: kb, only used when rule is ``twoClosest``
# -adv_oneDistance Unit: kb, only used when rule is ``oneClosest``
# -request_interval Time interval for two requests. Default is 300 seconds.
# -max_tries Maximum times trying to connect to GREAT web server.
# -version version of GREAT. The value should be "3.0.0", "2.0.2". Shorten version numbers
#          can also be used, such as using "3" or "3.0" is same as "3.0.0".
#
# == details
# Note it is not the standard GREAT API. This function directly send data to GREAT web server
# by HTTP POST.
#
# Following text is copied from GREAT web site ( http://bejerano-test.stanford.edu/great/public/html/index.php )
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
# == value
# A `GreatJob-class` class object which can be used to get results from GREAT server.
#
# == seealso
# `GreatJob-class`
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
submitGreatJob = function(gr, bg = NULL,
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
    version = "default"
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

    includeCuratedRegDoms = as.numeric(includeCuratedRegDoms[1])
    
    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")
    
    if(inherits(gr, "data.frame")) {
        gr = GRanges(seqnames = gr[[1]],
                     ranges = IRanges(start = gr[[2]],
                                       end = gr[[3]]))
    }
    gr = reduce(sort(gr))

    if(!bgChoice %in% c("wholeGenome", "data")) {
        stop("`bgChoice` should be one of 'wholeGenome' and 'data'.")
    }
    
    if(bgChoice == "data") {
        if(is.null(bg)) {
            stop("Since you set `bgChoice` to `data`, you need to provide values for `bg`.\n")
        }
    }
    if(!is.null(bg)) {
        if(inherits(bg, "data.frame")) {
            bg = GRanges(seqnames = bg[[1]],
                         ranges = IRanges(start = bg[[2]],
                                           end = bg[[3]]))
        }
        bg = reduce(sort(bg))
    }

    # check seqnames should have 'chr' prefix
    if(!all(grepl("^chr", seqnames(gr)))) {
        stop("Chromosome names (in `gr`) should have 'chr' prefix.\n")
    }
    if(!is.null(bg)) {
        if(!all(grepl("^chr", seqnames(bg)))) {
            stop("Chromosome names (in `bg`) should have 'chr' prefix.\n")
        }
    }

    # transform GRanges to data frame
    bed = as.data.frame(gr)
    bed_bg = NULL
    if(!is.null(bg)) {
        bed_bg = as.data.frame(bg)
    }

    
    # check request frequency
    time_interval = as.numeric(Sys.time()) - as.numeric(rGREAT_env$LAST_REQUEST_TIME)
    if(time_interval < request_interval) {
        message(qq("Don't make too frequent requests. The time break is @{request_interval}s.\nPlease wait for @{round(request_interval - time_interval)}s for the next request.\nThe time break can be set by `request_interval` argument.\n"))
        sleep(request_interval - time_interval)
    }

    #message("sending request to GREAT web server...")
    BASE_URL = BASE_URL_LIST[version]

    i_try = 0
    while(1) {
        
        error = try(response <- postForm(qq("@{BASE_URL}/greatWeb.php"),
                    "species"               = species,
                    "rule"                  = rule,
                    "span"                  = adv_span,
                    "upstream"              = adv_upstream,
                    "downstream"            = adv_downstream,
                    "twoDistance"           = adv_twoDistance,
                    "oneDistance"           = adv_oneDistance,
                    "includeCuratedRegDoms" = includeCuratedRegDoms,
                    "fgChoice"              = "data",
                    "fgData"                = qq("@{bed[[1]]}\t@{bed[[2]]}\t@{bed[[3]]}\t@{bed[[1]]}:@{bed[[2]]}-@{bed[[3]]}\n"),
                    "bgChoice"              = bgChoice,
                    "bgData"                = ifelse(bgChoice == "wholeGenome", "", qq("@{bed_bg[[1]]}\t@{bed_bg[[2]]}\t@{bed_bg[[3]]}\t@{bed_bg[[1]]}:@{bed_bg[[2]]}-@{bed_bg[[3]]}\n")),
                    "adv_upstream"          = adv_upstream,
                    "adv_downstream"        = adv_downstream,
                    "adv_span"              = adv_span,
                    "adv_twoDistance"       = adv_twoDistance,
                    "adv_oneDistance"       = adv_oneDistance
                    ))
        i_try = i_try + 1
        
        if(class(error) != "try-error") {
            break
        } else {
            message(error, appendLF = FALSE)
            if(i_try > max_tries) {
                stop(qq("max try: @{max_tries} reached. Stop with error.\n"))
            } else {
                message(qq("failed with the request, try after @{ceiling(request_interval/60)} min (@{i_try}th try)"))
                sleep(request_interval)
            }
        }
    }

    rGREAT_env$LAST_REQUEST_TIME = Sys.time()

    # parsing error
    if(any(grepl("encountered a user error", response))) {
        msg = gsub("^.*<blockquote>(.*?)<\\/blockquote>.*$", "\\1", response)
        msg = gsub("<.*?>", "", msg)
        msg = gsub(" +", " ", msg)
        msg = strwrap(msg)
        msg = paste(msg, collapse = "\n")
        stop(paste0("GREAT encountered a user error:\n", msg))
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
            "version"               = version)
    
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
    cat("Version:", param(object, "version"), "\n")
    #cat("Session ID:", id(object), "\n")
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

# == title
# Get enrichment tables from GREAT web server
#
# == param
# -job a `GreatJob-class` instance
# -ontology ontology names. Valid values are in `availableOntologies`. ``ontology`` is prior to 
#           ``category`` argument.
# -category Pre-defined ontology categories. One category can contain more than one ontologies. Valid values are in 
#            `availableCategories`
# -request_interval time interval for two requests. Default is 300 seconds.
# -max_tries maximum tries
#
# == details  
# The table contains statistics for the each term in each ontology catalogue.
#    
# Please note there is no FDR column in original tables. Users should 
# calculate by themselves by functions such as `stats::p.adjust`
# 
# == value
# The returned value is a list of data frames in which each one corresponds to 
# result for a single ontology. The structure of the data frames are same as 
# the tables available on GREAT website.
#
# == see also
# `availableOntologies`, `availableCategories`
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
setMethod(f = "getEnrichmentTables",
    signature = "GreatJob",
    definition = function(job, ontology = NULL, category = "GO",
    request_interval = 30, max_tries = 100) {
    
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
    
    res = lapply(ontology, function(onto) GREAT.read.json(job, qq("@{BASE_URL}/readJsFromFile.php?path=/scratch/great/tmp/results/@{jobid}.d/@{ONTOLOGY_KEYS[onto]}.js"), onto, 
        request_interval = request_interval, max_tries = max_tries))
    names(res) = ontology
    return(res)
})


# == title
# All available ontology names
#
# == param
# -job a `GreatJob-class` instance
# -category one or multiple categories. All available categories can be get by `availableCategories`
#
# == details
# The values of the supported ontologies sometime change. You should run the function to get the realtime
# values. The meaning of ontology returned is quite self-explained by the name.
#
# == value
# The returned values is a vector of ontologies.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
setMethod(f = "availableOntologies",
    signature = "GreatJob",
    definition = function(job, category = NULL) {
    
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
# Available ontology categories
#
# == param
# -job a `GreatJob-class` instance
#
# == details
# The values of the supported categories sometime change. You should run the function to get the realtime
# values. The meaning of categories returned is quite self-explained by the name.
#
# == value
# The returned value is a vector of categories.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
setMethod(f = "availableCategories",
    signature = "GreatJob",
    definition = function(job) {

    species = param(job, "species")
    CATEGORY = job@job_env$CATEGORY
    
    names(CATEGORY)
})

# download from `url` and save to `file`
download = function(url, file, request_interval = 30, max_tries = 100) {

    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")


    #message(qq("try to download from @{url}"))
    i_try = 0
    while(1) {
        error = try(download.file(url, destfile = file, quiet = TRUE))
            
        i_try = i_try + 1
            
        if(class(error) != "try-error") {
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

GREAT.read.json = function(job, url, onto, request_interval = 30, max_tries = 100) {
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
    
    if(!is.null(job@enrichment_tables[[onto]])) {
        res = job@enrichment_tables[[onto]]
        return(res)
    }

    download(url, file = f1, request_interval = request_interval, max_tries = max_tries)

    # if something wrong, the returned js has length of 0
    if(file.info(f1)$size == 0) {
        stop("Retrieved 0 byte from remote server. Probably your job on GREAT server expires\nand you need to submit it again.\n")
    }
    
    # just in case downloading json file is interrupted
    error = try(json <- fromJSON(file = f1))
    if(class(error) == "try-error") {
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
    
    colnames(res) = c("ID", "name", "Binom_Genome_Fraction", "Binom_Expected", "Binom_Observed_Region_Hits", "Binom_Fold_Enrichment",
                      "Binom_Region_Set_Coverage", "Binom_Raw_PValue", "Hyper_Total_Genes", "Hyper_Expected",
                      "Hyper_Observed_Gene_Hits", "Hyper_Fold_Enrichment", "Hyper_Gene_Set_Coverage",
                      "Hyper_Term_Gene_Coverage", "Hyper_Raw_PValue")
    job@enrichment_tables[[onto]] = res
    return(res)
}

# transform the original file into the standard data frame
parseRegionGeneAssociationFile = function(f1) {

    error = try(data <- read.table(f1, sep = "\t", stringsAsFactors = FALSE))
    if(class(error) == "try-error") {
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
            r = gregexpr("([a-zA-Z0-9\\-_.]+) \\(([+-]\\d+)\\)", data[i, 2], perl = TRUE)[[1]]
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
# Plot region-gene association figures
#
# == param
# -job a `GreatJob-class` instance
# -type type of plots, should be in ``1, 2, 3``. See details section for explanation
# -ontology ontology name
# -termID term id which corresponds to the selected ontology
# -request_interval time interval for two requests. Default is 300 seconds.
# -max_tries maximum tries
#
# == details
# Generated figures are:  
#
# - association between regions and genes
# - distribution of distance to TSS
# - distribution of absolute distance to TSS
#
#     
# If ``ontology`` and ``termID`` are set, only regions and genes corresponding to 
# selected ontology term will be used. Valid value for ``ontology`` is in 
# `availableOntologies` and valid value for ``termID`` is from 'id' column 
# in the table which is returned by `getEnrichmentTables`.  
#
# == value
# a `GenomicRanges::GRanges` object. Columns in metadata are:  
#
# -gene genes that are associated with corresponding regions
# -distTSS distance from the regions to TSS of the associated gene
#
# The returned values corresponds to whole input regions or only regions in specified ontology term, 
# depending on user's setting. 
#
# If there is no gene associated with the region, corresponding ``gene`` and ``distTSS``
# columns will be ``NA``.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
setMethod(f = "plotRegionGeneAssociationGraphs",
    signature = "GreatJob",
    definition = function(job, type = 1:3, ontology = NULL, 
    termID = NULL, request_interval = 30, max_tries = 100) {

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
        stop("You should set both of `ontology` and `termID` or neither of them.\n")
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
            url = qq("@{BASE_URL}/downloadAssociations.php?termId=@{termID}&ontoName=@{ONTOLOGY_KEYS[ontology]}&sessionName=@{jobid}&species=@{species}&foreName=user-provided%20data&backName=&table=region")
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
        url = qq("@{BASE_URL}/downloadAssociations.php?sessionName=@{jobid}&species=@{species}&foreName=user-provided%20data&backName=&table=region")
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

    # make plots
    if(1 %in% type) {
        if(using_term) {
            tb = table(table(df_term$gene))
            vt = numeric(11)
            vt[as.numeric(names(tb))] = tb
            vt[is.na(vt)] = 0
            v = c(vt[1:10], sum(vt[10:length(vt)]))
            names(v) = c(as.character(1:10), ">10")
            p = v/sum(v)
            pos = barplot(p, col = "black", xlab = "Number of associated regions per gene", ylab = "This term's genes", ylim = c(0, max(p)*1.5), main = qq("Number of associated regions per gene\n@{ontology}\n@{termID}"))
            text(pos[, 1], p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        } else {
            tb = table(table(paste(df_all$chr, df_all$start, df_all$end, sep = ",")))
            v = c(nrow(df_all_NA), tb["1"], tb["2"], sum(tb[as.numeric(names(tb)) > 2]))
            names(v) = c("0", "1", "2", "> 3")
            p = v/sum(v)
            pos = barplot(p, col = c("red", "grey", "grey", "grey"), xlab = "Number of associated genes per region", ylab = "Genomic regions", ylim = c(0, max(p)*1.5), main = "Number of associated genes per region")
            text(pos[, 1], p + 0.01, v, adj = c(0.5, 0), col = c("red", "black", "black", "black"), cex = 0.8)
            legend("topright", pch = 15, col = c("grey", "red"), legend = c("Genomic regions associated with one or more genes", "Genomic regions not associated with any genes"), cex = 0.8)
        }
    }
    if(2 %in% type) {
        v = cbind(
            c("<-500"       = sum(df_all$distTSS <= -500000),
              "-500 to -50" = sum(df_all$distTSS > -500000 & df_all$distTSS <= -50000),
              "-50 to -5"   = sum(df_all$distTSS > -50000  & df_all$distTSS <= -5000),
              "-5 to 0"     = sum(df_all$distTSS > -5000   & df_all$distTSS <= 0),
              "0 to 5"      = sum(df_all$distTSS > 0       & df_all$distTSS <= 5000),
              "5 to 50"     = sum(df_all$distTSS > 5000    & df_all$distTSS <= 50000),
              "50 to 500"   = sum(df_all$distTSS > 50000   & df_all$distTSS <= 500000),
              "> 500"       = sum(df_all$distTSS > 500000)))
        p = v/sum(v)
        if(using_term) {
            v = cbind( 
            c("<-500"       = sum(df_term$distTSS <= -500000),
              "-500 to -50" = sum(df_term$distTSS > -500000 & df_term$distTSS <= -50000),
              "-50 to -5"   = sum(df_term$distTSS > -50000  & df_term$distTSS <= -5000),
              "-5 to 0"     = sum(df_term$distTSS > -5000   & df_term$distTSS <= 0),
              "0 to 5"      = sum(df_term$distTSS > 0       & df_term$distTSS <= 5000),
              "5 to 50"     = sum(df_term$distTSS > 5000    & df_term$distTSS <= 50000),
              "50 to 500"   = sum(df_term$distTSS > 50000   & df_term$distTSS <= 500000),
              "> 500"       = sum(df_term$distTSS > 500000)), v)
            p = v/sum(v)
        }
        
        rownames(p) = NULL
        if(using_term) {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "Distance to TSS (kb)", ylab = "Region-gene associations", ylim = c(0, max(p)*1.5), main = qq("Binned by orientation and distance to TSS\n@{ontology}\n@{termID}"), axes = FALSE)
        } else {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "Distance to TSS (kb)", ylab = "Region-gene associations", ylim = c(0, max(p)*1.5), main = qq("Binned by orientation and distance to TSS"), axes = FALSE)
        }
        text(t(pos), p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        axis(side = 2)
        op = par("xpd")
        par(xpd = NA)
        text(colMeans(pos), -0.02,rownames(v), srt = 45, adj = c(1, 0.5))
        par(xpd = op)
        x1 = (mean(pos[, 4]) + mean(pos[, 5]))/2
        x2 = (mean(pos[, 5]) + mean(pos[, 6]))/2
        y = max(p)*1.2
        lines( c(x1, x1), c(0, y))
        arrows(x1, y, x2, y, angle = 15, length = 0.1, code = 2)
        text(mean(pos[, 5]), y+0.01, "TSS", adj = c(0.5, 0), cex = 0.8)
        if(using_term) {
            legend("topright", pch = 15, col = c("green", "blue"), legend = c("This term's region-gene associations", "Set-wide region-gene associations"), cex = 0.8)
        }
    }
    if(3 %in% type) {
        v = cbind(
            c("0 to 5"      = sum(abs(df_all$distTSS) > 0       & abs(df_all$distTSS) <= 5000),
              "5 to 50"     = sum(abs(df_all$distTSS) > 5000    & abs(df_all$distTSS) <= 50000),
              "50 to 500"   = sum(abs(df_all$distTSS) > 50000   & abs(df_all$distTSS) <= 500000),
              "> 500"       = sum(abs(df_all$distTSS) > 500000)))
        p = v/sum(v)
        if(using_term) {
            v = cbind( 
            c("0 to 5"      = sum(abs(df_term$distTSS) > 0       & abs(df_term$distTSS) <= 5000),
              "5 to 50"     = sum(abs(df_term$distTSS) > 5000    & abs(df_term$distTSS) <= 50000),
              "50 to 500"   = sum(abs(df_term$distTSS) > 50000   & abs(df_term$distTSS) <= 500000),
              "> 500"       = sum(abs(df_term$distTSS) > 500000)), v)
            p = v/sum(v)
        }
        
        rownames(p) = NULL
        if(using_term) {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "Absolute distance to TSS (kb)", ylab = "Region-gene associations", ylim = c(0, max(p)*1.5), main = qq("Binned by absolute distance to TSS\n@{ontology}\n@{termID}"), axes = FALSE)
        } else {
            pos = barplot(t(p), beside = TRUE, col = {if(using_term) c("green", "blue") else "blue"}, xlab = "Absolute distance to TSS (kb)", ylab = "Region-gene associations", ylim = c(0, max(p)*1.5), main = qq("Binned by absolute distance to TSS"), axes = FALSE)
        }
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
            legend("topright", pch = 15, col = c("green", "blue"), legend = c("This term's region-gene associations", "Set-wide region-gene associations"))
        }
    }
    
    if(using_term) {
        df = rbind(df_term, df_term_NA)
    } else {
        df = rbind(df_all, df_all_NA)
    }

    gr = GRanges(seqnames = factor(df[[1]], levels = sort_chr(unique(df[[1]]))),
                 ranges = IRanges(start = df[[2]],
                                   end = df[[3]]),
                 gene = df[[4]],
                 distTSS = df[[5]])
    gr = sort(gr)
    return(invisible(gr))
})

# just want to print something while waiting
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
    if(file.info(file)$size < 10000) {

        # if session expires, there is a error page
        file_text = readLines(file)
        if(any(grepl("GREAT is unable to locate a session", file_text))) {
            stop("Your job on GREAT server expires and you need to submit it again.\n")
        }
        
        # wrong termID returns a table only contain a newline
        if(all(grepl("^\\s*$", file_text))) {
            stop("Empty data, probably your 'termID' is invalid.\n")
        }
    }


}
