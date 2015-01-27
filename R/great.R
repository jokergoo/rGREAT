
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
#
# == details
# Note it is not the standard GREAT API. This function directly send data to GREAT web server
# by HTTP POST.
#
# Following text is copied from GREAT web site ( http://bejerano-test.stanford.edu/great/public/html/index.php )
#
# Explanation of ``rule`` and settings with names started with 'adv_':
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
# A ``GREAT_Job`` class object which can be used to get results from GREAT server.
#
# == seealso
# `GREAT_Job-class`
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
submitGREATJob = function(gr, bg = NULL,
    species               = c("hg19", "hg18", "mm9", "danRer7"),
    includeCuratedRegDoms = TRUE,
    bgChoice              = c("wholeGenome", "data"),
    rule                  = c("basalPlusExt", "twoClosest", "oneClosest"),
    adv_upstream          = 5.0,
    adv_downstream        = 1.0,
    adv_span              = 1000.0,
    adv_twoDistance       = 1000.0,
    adv_oneDistance       = 1000.0,
	request_interval = 300,
	max_tries = 10
    ) {
    
    species  = match.arg(species)[1]
    rule     = match.arg(rule)[1]
    bgChoice = match.arg(bgChoice)[1]
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
        message(qq("Don't make too frequent requests. The time break is @{request_interval}s.\nPlease wait for @{round(request_interval - time_interval)}s for the next request.\nThe time break can be set by `GREAT.options()`.\n"))
        sleep(request_interval - time_interval)
    }
    
    #message("sending request to GREAT web server...")
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
                    "bgData"                = ifelse(bgChoice == "wholeGenome", "", qq("@{bed_bg[[1]]}\t@{bed_bg[[2]]}\t@{bed_bg[[3]]}\n")),
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
   
    if(any(grepl("encountered a user error", response))) {
        stop("GREAT encountered a user error, check your input (especially `species`).\n")
    }
    
    td = tempdir()
    dir.create(td, showWarnings = FALSE)
    
    jobid = gsub("^.*var _sessionName = \"(.*?)\";.*$", "\\1",  response)
    jobid = as.vector(jobid)
      
    job = GREAT_Job(id = jobid,
        parameters = list(
            submit_time = Sys.time(),
            species = species,
            bgChoice = bgChoice,
            includeCuratedRegDoms = includeCuratedRegDoms,
            rule = rule,
            adv_upstream = adv_upstream,
            adv_downstream = adv_downstream,
            adv_span = adv_span,
            adv_twoDistance = adv_twoDistance,
            adv_oneDistance = adv_oneDistance,
            tempdir = td))
    
    return(job)
}


GREAT_Job$methods(getEnrichmentTables = function(ontology = NULL, category = c("GO", "Pathway_Data"),
	request_interval = 30, max_tries = 100) {
    
    jobid = .self$get_id()
    species = .self$get_param("species")
    
    if(is.null(ontology)) {
        if(is.null(category)) {
            stop("`ontology` and `category` can not be both NULL.\n")
        }
        if(length(setdiff(category, .self$availableCategories()))) {
            stop("Only categories in `job$availableCategories()` are allowed.\n")
        }
        ontology = .self$availableOntologies(category = category)
    }
    ontology = unique(ontology)
    
    if(length(setdiff(ontology, .self$availableOntologies()))) {
        stop("Only ontologies in `job$availableOntologies()` are allowed.\n")
    }
    
    res = lapply(ontology, function(onto) GREAT.read.json(job, qq(URL_TEMPLATE[onto]), onto, 
		request_interval = request_interval, max_tries = max_tries))
    names(res) = ontology
    return(res)
})


GREAT_Job$methods(availableOntologies = function(category = NULL) {
    
    species = .self$get_param("species")
    
    if(is.null(category)) {
        onto = unlist(CATEGORY[[species]])
        names(onto) = NULL
        return(onto)
    } else {
        if(length(setdiff(category, .self$availableCategories())) > 0) {
            stop("Value of `category` is invalid. Please use job$availableCategories() to find supported categories.\n")
        }
        onto = unlist(CATEGORY[[species]][category])
        names(onto) = NULL
        return(onto)
    }
})


GREAT_Job$methods(availableCategories = function() {

    species = .self$get_param("species")
    
    names(CATEGORY[[species]])
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
    jobid = job$get_id()
    TEMP_DIR = job$get_param("tempdir")

    op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")

    
    f1 = qq("@{TEMP_DIR}/@{jobid}_@{onto}.js")
    
    if(!is.null(job$enrichment_tables[[onto]])) {
        res = job$enrichment_tables[[onto]]
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
    job$enrichment_tables[[onto]] = res
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


GREAT_Job$methods(plotRegionGeneAssociationGraphs = function(type = 1:3, ontology = NULL, 
	termID = NULL, request_interval = 30, max_tries = 100) {

    jobid = .self$get_id()
    species = .self$get_param("species")
    TEMP_DIR = .self$get_param("tempdir")

    opqq = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(opqq))
    qq.options(code.pattern = "@\\{CODE\\}")

    # make plot and return value for a single term
    using_term = FALSE
    
    if(!is.null(ontology) && !is.null(termID)) {
        using_term = TRUE
        
        ontology = ontology[1]
        termID = termID[1]

        if(! ontology %in% .self$availableOntologies()) {
			stop("Value of `ontology` should be in `job$availableOntologies()`\n")
		}
    }

    if(sum(c(is.null(ontology), is.null(termID))) == 1) {
        stop("You should set both of `ontology` and `termID` or neither of them.\n")
    }
	
	if(using_term) {
		# check whether termID is in ontology if ontology table is already downloaded
		if(!is.null(job$enrichment_tables[[ontology]])) {
			if(! termID %in% job$enrichment_tables[[ontology]]$ID) {
				stop(qq("Cannot find '@{termID}' in enrichment table of '@{ontology}'"))
			}
		}
	}
    
    if(using_term) {

        # prepare file names for local table
        f_term = qq("@{jobid}-@{ontology}-@{termID}-@{species}-region.txt")
        f_term = gsub('[\\/:*?"<>|]', "_", f_term)
        f_term = qq("@{TEMP_DIR}/@{f_term}")
        
        if(!is.null(job$association_tables[[qq("@{ontology}-@{termID}")]])) {
            df_term = job$association_tables[[qq("@{ontology}-@{termID}")]]
        } else {
            url = qq("@{BASE_URL}/downloadAssociations.php?termId=@{termID}&ontoName=@{ONTOLOGY_KEYS[ontology]}&sessionName=@{jobid}&species=@{species}&foreName=user-provided%20data&backName=&table=region")
            download(url, file = f_term, request_interval = request_interval, max_tries = max_tries)
            check_asso_file(f_term)
            df_term = parseRegionGeneAssociationFile(f_term)
            job$association_tables[[qq("@{ontology}-@{termID}")]] = df_term
            file.remove(f_term)
        }
    }
    
    f_all = qq("@{TEMP_DIR}/@{jobid}-@{species}-all-region.txt")
        
    # download
    if(!is.null(job$association_tables[["all"]])) {
        df_all = job$association_tables[["all"]]
    } else {
        url = qq("@{BASE_URL}/downloadAssociations.php?sessionName=@{jobid}&species=@{species}&foreName=user-provided%20data&backName=&table=region")
        download(url, file = f_all, request_interval = request_interval, max_tries = max_tries)
        check_asso_file(f_all)
        df_all = parseRegionGeneAssociationFile(f_all)
        job$association_tables[["all"]] = df_all
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
            p = v/(nrow(df_term) + nrow(df_term_NA))
            pos = barplot(p, col = "black", xlab = "Number of associated regions per gene", ylab = "This term's genes", ylim = c(0, max(p)*1.5), main = qq("Number of associated regions per gene\n@{ontology}\n@{termID}"))
            text(pos[, 1], p + 0.01, v, adj = c(0.5, 0), cex = 0.8)
        } else {
            tb = table(table(paste(df_all$chr, df_all$start, df_all$end, sep = ",")))
            v = c(nrow(df_all_NA), tb["1"], tb["2"], sum(tb[as.numeric(names(tb)) > 2]))
            names(v) = c("0", "1", "2", "> 3")
            p = v/(nrow(df_all) + nrow(df_all_NA))
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
        p = v/(nrow(df_all) + nrow(df_all_NA))
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
            p = v/c(rep(nrow(df_term) + nrow(df_term_NA), nrow(v)), rep(nrow(df_all) + nrow(df_all_NA), nrow(v)))
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
        p = v/(nrow(df_all) + nrow(df_all_NA))
        if(using_term) {
            v = cbind( 
            c("0 to 5"      = sum(abs(df_term$distTSS) > 0       & abs(df_term$distTSS) <= 5000),
              "5 to 50"     = sum(abs(df_term$distTSS) > 5000    & abs(df_term$distTSS) <= 50000),
              "50 to 500"   = sum(abs(df_term$distTSS) > 50000   & abs(df_term$distTSS) <= 500000),
              "> 500"       = sum(abs(df_term$distTSS) > 500000)), v)
            p = v/c(rep(nrow(df_term) + nrow(df_term_NA), nrow(v)), rep(nrow(df_all) + nrow(df_all_NA), nrow(v)))
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