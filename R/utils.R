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


get_table_from_orgdb = function(pattern, orgdb) {
	check_pkg(orgdb)
	all_tb = ls(envir = getNamespace(orgdb))
	i_tb = which(grepl(pattern, all_tb))
	if(length(i_tb)) {
		getFromNamespace(all_tb[i_tb], orgdb)
	} else {
		NULL
	}	
}

# == title
# Generate random regions
#
# == param
# -genome UCSC genome version, e.g. "hg19".
# -nr Number of regions.
# -seqlengths Alternatively, you can also specify a named vector of seqlengths (chromosome lengths).
# -width_fun A function which defines the distribution of region widths.
#
# == details
# The number of regions per chromosome is proportional to the chromsome length.
#
# == example
# gr = randomRegions(genome = "hg19")
# quantile(width(gr))
randomRegions = function (genome = NULL, nr = 1000, seqlengths = NULL,
	width_fun = function(n) runif(n, min = 1000, max = 10000)) {

	if(is.null(genome) && is.null(seqlengths)) {
		stop_wrap("Either `genome` or `seqlengths` should be specified.")
	}

	if(!is.null(genome)) {
	    chr.len = getChromInfoFromUCSC(genome)
	    chromosome = names(chr.len)
	} else {
		chr.len = seqlengths
		chromosome = names(chr.len)
	}

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

    gr = GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df[, 3]))
    seqlevels(gr) = unique(chromosome)
    # if(!is.null(genome)) genome(gr) = genome

    return(gr)
}

guess_gene_id_type = function(id) {
    l = grepl("^\\d+$", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENTREZ")
    }
    l = grepl("^ENS.*G", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBL")
    }
    l = grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        return("REFSEQ")
    }
    return(NULL)
}


# == title
# Read gmt gene sets file
#
# == param
# -x The file name of a .gmt file.
# -from Gene ID type in the original gmt file. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.
# -to Gene ID type that you want to convert to. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.
# -orgdb The name of an OrgDb database.
#
# == value
# A named list of vectors.
#
# == example
# read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"))
read_gmt = function(x, from = NULL, to = NULL, orgdb = NULL) {
	ln = readLines(x)
	lt = strsplit(ln, "\t")

	nm = sapply(lt, function(x) x[1])
	v = lapply(lt, function(x) x[-(1:2)])
	names(v) = nm

	if(!is.null(from) && !is.null(to) && !is.null(orgdb)) {

		if(!from %in% c("ENTREZ", "SYMBOL", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`from` can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.")
		}
		if(!to %in% c("ENTREZ", "SYMBOL", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`to` can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.")
		}

		if(from == "ENTREZ") {
			map = get_table_from_orgdb(paste0(to, "$"), orgdb)
			map = unlist(as.list(map))
			v = lapply(v, function(x) {
				x2 = map[x]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		} else if(to == "ENTREZ") {
			map = get_table_from_orgdb(paste0(from, "2EG$"), orgdb)
			map = unlist(as.list(map))
			v = lapply(v, function(x) {
				x2 = map[x]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		} else {
			map1 = get_table_from_orgdb(paste0(from, "2EG$"), orgdb)
			map1 = unlist(as.list(map1))

			map2 = get_table_from_orgdb(paste0(to, "$"), orgdb)
			map2 = unlist(as.list(map2))
			v = lapply(v, function(x) {
				x2 = map2[map1[x]]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		}
	}
	v = v[sapply(v, length) > 0]
	v
}


get_url = function(url) {
	nm = basename(url)

	f = paste0(tempdir(), "/", nm)
	if(!file.exists(f)) {
		e = try(suppressWarnings(download.file(url, destfile = f, quiet = TRUE)), silent = TRUE)
		if(inherits(e, "try-error")) {
			stop_wrap(qq("Cannot download file from @{url}."))
		}
	}
	f
}

# chr_len_db = tapply(BIOC_ANNO_PKGS$txdb, BIOC_ANNO_PKGS$genome_version_in_txdb, function(x) {
# 	x = x[1]
# 	check_pkg(x, bioc = TRUE)
# 	txdb = getFromNamespace(x, ns = x)
# 	seqlengths(txdb)
# })

# chr_len_db = lapply(chr_len_db, function(x) {
# 	x = x[nchar(names(x)) < 10]
# 	if(grepl("$chr", names(x)[1])) {
# 		names(x) = paste0("chr", names(x))
# 	}
# 	x
# })
