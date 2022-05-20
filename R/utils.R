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
# Get RefSeq genes from UCSC
#
# == param
# -genome UCSC genome, such as "hg19".
# -subset Subset of RefSeq genes. See https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=refSeqComposite .
#
# == value
# A `GenomicRanges` object.
getRefSeqGenesFromUCSC = function(genome, subset = c("RefSeqSelect", "RefSeqCurated", "RefSeq")) {
	subset = match.arg(subset)[1]
	url = qq("https://hgdownload.cse.ucsc.edu/goldenPath/@{genome}/database/ncbi@{subset}.txt.gz")

	refseq = paste0(tempdir(), qq("/@{genome}_ncbi@{subset}.txt.gz"))
	if(!file.exists(refseq)) {
		e = try(suppressWarnings(download.file(url, destfile = refseq, quiet = TRUE)), silent = TRUE)
		if(inherits(e, "try-error")) {
			stop_wrap(qq("It seems UCSC does not provide 'ncbi@{genome}.txt.gz' for @{genome} or the internet connection was interrupted. If possible, download gap file directly from @{url}."))
		}
	}
	
	tb = read.table(gzfile(refseq), sep = "\t", stringsAsFactors = FALSE)[, c(3, 5, 6, 4, 2)]
	tb[, 5] = gsub("\\.\\d+", "", tb[, 5])
	tb[, 2] = tb[, 2] + 1

	gene_id_type = "REFSEQ"

	i = which(BIOC_ANNO_PKGS$genome_version_in_txdb == genome)
	if(length(i)) {
		orgdb = BIOC_ANNO_PKGS$orgdb[i][1]
		map = get_table_from_org_db("REFSEQ2EG$", orgdb)
		map = unlist(as.list(map))
		tb[, 5] = map[tb[, 5]]
		tb = tb[!is.na(tb[, 5]), ,drop = FALSE]
		gene_id_type = "ENTREZ"
	}
	
	gr = GRanges(seqnames = tb[, 1], ranges = IRanges(tb[, 2], tb[, 3]), strand = tb[, 4], gene_id = tb[, 5])
	gr = unique(gr)

	m = metadata(gr)
	attr(m, "genome") = genome
	attr(m, "gene_id_type") = gene_id_type
	metadata(gr) = m
	gr
}


# == title
# Get built-in TSS from GREAT
#
# == param
# -genome Only support "hg19", "hg38", "mm10", "mm9". Files are downloaded from https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655445/Genes .
#
# == value
# A `GenomicRanges::GRanges` object.
getGREATDefaultTSS = function(genome) {
	if(!genome %in% c("hg19", "hg38", "mm9", "mm10")) {
		stop_wrap("`genome` can only take value in 'hg19', 'hg38', 'mm9', 'mm10'.")
	}

	tb = read.table(system.file("extdata", paste0("GREATv4.genes.", genome, ".tsv.gz"), package = "rGREAT"))

	gr = GRanges(seqnames = tb[, 2], ranges = IRanges(tb[, 3], tb[, 3]), strand = tb[, 4], gene_id = tb[, 5])

	m = metadata(gr)
	attr(m, "genome") = genome
	attr(m, "gene_id_type") = "SYMBOL"
	metadata(gr) = m
	gr
}

# == title
# Get Gencode genes
#
# == param
# -version Gencode version, e.g. v19 for human, vM21 for mouse.
#
# == details
# Only the protein coding genes.
#
# == value
# A `GenomicRanges::GRanges` object.
getGencodeGenes = function(version) {
	lt = readRDS(system.file("extdata", "gencode_gene.rds", package = "rGREAT"))

	if(!qq("gencode.@{version}.annotation.gtf.gz") %in% names(lt)) {
		if(qq("gencode.v@{version}.annotation.gtf.gz") %in% names(lt)) {
			version = paste0("v", version)
		} else {
			stop_wrap(qq("Cannot find gencode version @{version}."))
		}
	}

	gr = lt[[qq("gencode.@{version}.annotation.gtf.gz")]]

	if(version %in% paste0("v", 4:19)) {
		genome = "hg19"
	} else if(version %in% paste0("v", 20:40)) {
		genome = "hg38"
	} else if(version %in% paste0("vM", 1)) {
		genome = "mm9"
	} else if(version %in% paste0("vM", 2:25)) {
		genome = "mm10"
	} else if(version %in% paste0("vM", 26:29)) {
		genome = "mm11"
	}

	m = metadata(gr)
	attr(m, "genome") = genome
	attr(m, "gene_id_type") = "ENSEMBL"
	metadata(gr) = m
	gr
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
