check_pkg = function(pkg, bioc = FALSE) {
	if(requireNamespace(pkg, quietly = TRUE)) {
		return(NULL)
	} else {

		if(!interactive()) {
			if(bioc) {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}

		if(bioc) {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from Bioconductor? [y|n] "))
		} else {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from CRAN? [y|n] "))
		}

		if(bioc) {
			if(tolower(answer) %in% c("y", "yes")) {
				if(!requireNamespace("BiocManager", quietly = TRUE)) {
					install.packages("BiocManager")
				}
				BiocManager::install(pkg)
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			}
		} else {
			if(tolower(answer) %in% c("y", "yes")) {
				install.packages(pkg)
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}
	}
}


stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

warning_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    warning(x, call. = FALSE)
}

message_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}


get_table_from_org_db = function(pattern, org_db) {
	check_pkg(org_db)
	all_tb = ls(envir = getNamespace(org_db))
	i_tb = which(grepl(pattern, all_tb))
	if(length(i_tb)) {
		getFromNamespace(all_tb[i_tb], org_db)
	} else {
		NULL
	}	
}

# == title
# Get gap regions from UCSC
#
# == param
# -genome UCSC genome, such as "hg19".
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == example
# getGapFromUCSC("hg19")
getGapFromUCSC = function(genome) {

	url = paste0("https://hgdownload.cse.ucsc.edu/goldenPath/", genome, "/database/gap.txt.gz")


	gap = paste0(tempdir(), "/", genome, "_gap.txt.gz")
	if(!file.exists(gap)) {
		e = try(suppressWarnings(download.file(url, destfile = gap, quiet = TRUE)), silent = TRUE)
		if(inherits(e, "try-error")) {
			stop_wrap(qq("It seems UCSC does not provide 'gap.txt.gz' for @{genome} or the internet connection was interrupted. If possible, download gap file directly from @{url}."))
		}
	}
	
	tb = read.table(gzfile(gap), sep = "\t", stringsAsFactors = FALSE)[, 2:4]
	tb[, 2] = tb[, 2] + 1
	GRanges(seqnames = tb[, 1], ranges = IRanges(tb[, 2], tb[, 3]))
}


# == title
# Generate random regions
#
# == param
# -nr Number of regions.
# -genome UCSC genome version, e.g. "hg19".
# -width_fun A function which defines the distribution of region widths.
#
# == details
# The number of regions per chromosome is proportional to the chromsome length.
#
# == example
# gr = randomRegions(1000)
# quantile(width(gr))
randomRegions = function (nr = 1000, genome = NULL,
	width_fun = function(n) runif(n, min = 1000, max = 10000)) {

    cyto = read.cytoband(species = genome)
    chr.len = cyto$chr.len
    chromosome = cyto$chromosome
    dl = lapply(seq_along(chr.len), function(i) {
        k = round(nr * chr.len[i]/sum(chr.len))
        breaks = sort(sample(chr.len[i], k))
        
        rr = data.frame(chr = rep(chromosome[i], length(breaks)),
        	       start = breaks,
        	       end = breaks + round(width_fun(length(breaks))))

        rr = rr[rr[, 3] <= chr.len[i], , drop = FALSE]
        rr
    })
    df = NULL
    for (i in seq_along(dl)) {
        df = rbind(df, dl[[i]])
    }
    
    return(GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df[, 3])))
}
