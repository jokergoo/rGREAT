

if(identical(topenv(), .GlobalEnv)) {
	BIOC_ANNO_PKGS = read.table("~/project/development/rGREAT/inst/extdata/bioc_anno_pkgs.csv", header = TRUE, sep = ";")
	CHR_LEN_DB = readRDS("~/project/development/rGREAT/inst/extdata/CHR_LEN_DB.rds")
} else {
	BIOC_ANNO_PKGS = read.table(system.file("extdata", "bioc_anno_pkgs.csv", package = "rGREAT"), header = TRUE, sep = ";")
	CHR_LEN_DB = readRDS(system.file("extdata", "CHR_LEN_DB.rds", package = "rGREAT"))
}

rGREAT_env$BIOMART = NULL

rGREAT_env$extended_tss = list()

detect_txdb = function(x) {
	i = which(BIOC_ANNO_PKGS$txdb %in% x)
	if(length(i)) {
		return(i)
	}

	x = tolower(as.character(x))

	# hg19
	i = which(tolower(BIOC_ANNO_PKGS$genome_version_in_txdb) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# human
	i = which(tolower(BIOC_ANNO_PKGS$species_name) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# Homo sapiens
	i = which(tolower(BIOC_ANNO_PKGS$species_latin) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# 9606
	i = which(as.character(BIOC_ANNO_PKGS$taxon_id) %in% x)
	if(length(i)) {
		return(i[1])
	}

	return(i)
}

# == title
# Extend TSS
#
# == param
# -txdb Name of "TxDb.*" packages from Bioconductor. All supported TxDb packages are in ``rGREAT:::BIOC_ANNO_PKGS$txdb``. Note short genome version can also be used
#     here such as "hg19" or "hg19.knownGene".
# -verbose Whether to print messages.
# -... All pass to `extendTSS`.
#
# == value
# A `GenomicRanges::GRanges` object with one meta column 'gene_id'.
#
# == example
# if(FALSE) {
# extendTSSFromTxDb("TxDb.Hsapiens.UCSC.hg19.knownGene")
# extendTSSFromTxDb("hg19")
# }
extendTSSFromTxDb = function(txdb, verbose = great_opt$verbose, ...) {

	txdb_pkg = txdb
	if(!inherits(txdb_pkg, "character")) {
		txdb_pkg = txdb_pkg$packageName
	}

	i = detect_txdb(txdb_pkg)
	if(length(i) == 0) {
		stop_wrap(qq("TxDb package '@{txdb_pkg}' is not supported."))
	}

	txdb_pkg = BIOC_ANNO_PKGS$txdb[i]

	hash = digest::digest(list(txdb_pkg, ...))

	if(!is.null(rGREAT_env$extended_tss[[hash]])) {
		if(verbose) {
			message("* extended_tss is already cached, directly use it.")
		}
		return(rGREAT_env$extended_tss[[hash]])
	}
	
	if(verbose) {
		message(qq("* check whether TxDb package '@{txdb_pkg}' is installed."))
	}

	check_pkg(txdb_pkg, bioc = TRUE)

	if(verbose) {
		message(qq("* gene ID type in the extended TSS is '@{BIOC_ANNO_PKGS$gene_id_in_txdb[i]}'."))
	}

	suppressMessages(gene <- genes( getFromNamespace(txdb_pkg, txdb_pkg) ))

	genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i]
	cl = CHR_LEN_DB[[genome]]
	chromosomes = names(cl)

	if(verbose) {
		message_str = qq("* restrict chromosomes to '@{paste(chromosomes, collapse=', ')}'.")
		message_str = strwrap(message_str, width = 90)
		if(length(message_str) > 1) {
			message_str[2:length(message_str)]= paste0("    ", message_str[2:length(message_str)])
		}
		message(paste(message_str, collapse = "\n"))
	}

	gene = gene[seqnames(gene) %in% chromosomes]

	orgdb = BIOC_ANNO_PKGS$orgdb[i]

	check_pkg(orgdb, bioc = TRUE)

	all_tb = ls(envir = getNamespace(orgdb))
	i_tb = which(grepl("ENSEMBL2EG", all_tb))

	gene_id_type = BIOC_ANNO_PKGS$gene_id_in_txdb[i]
	if(BIOC_ANNO_PKGS$gene_id_in_txdb[i] == "Ensembl gene ID" & length(i_tb)) {
		if(verbose) {
			message("* convert gene ID from 'Ensembl gene ID' to 'Entrez gene ID'.")
		}
		map = as.list(get_table_from_orgdb("ENSEMBL2EG", orgdb))
		map = sapply(map, function(x) x[1])
		names(map) = gsub("\\.\\d+$", "", names(map))

		gene_id = names(gene)
		new_gene_id = map[ gsub("\\.\\d+$", "", gene_id) ]

		l = !is.na(new_gene_id)

		gene = gene[l, ]
		gene$gene_id = new_gene_id[l]
		names(gene) = new_gene_id[l]

		gene_id_type = "Entrez Gene ID"
	}

	i_tb2 = which(grepl("GENETYPE", all_tb))
	if(length(i_tb2)) {
		gene_type_tb = getFromNamespace(all_tb[i_tb2], orgdb)
		gene_type = unlist(as.list(gene_type_tb))
		ind = gene_type[names(gene)] == "protein-coding"
		ind[is.na(ind)] = FALSE
		gene = gene[ind]

		if(verbose) {
			message(qq("* @{length(gene)}/@{length(ind)} protein-coding genes left."))
		}
	} else {
		if(verbose) {
			message(qq("* @{length(gene)} genes left."))
		}
	}

	if(verbose) {
		message(qq("* update seqinfo to the selected chromosomes."))
	}
	seqlevels(gene) = chromosomes
	info = seqinfo(gene)
	info = info[chromosomes, ]
	seqinfo(gene) = info

	extended_tss = extendTSS(gene, seqlengths(gene), verbose = verbose, 
		.attr = list(genome = genome, gene_id_type = gene_id_type, txdb = txdb_pkg, orgdb = orgdb), ...)

	rGREAT_env$extended_tss[[hash]] = extended_tss

	return(extended_tss)
}


detect_orgdb = function(x) {
	i = which(BIOC_ANNO_PKGS$orgdb %in% x)
	if(length(i)) {
		return(i[1])
	}

	x = tolower(as.character(x))

	# hg19
	i = which(tolower(BIOC_ANNO_PKGS$genome_version_in_txdb) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# human
	i = which(tolower(BIOC_ANNO_PKGS$species_name) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# Homo sapiens
	i = which(tolower(BIOC_ANNO_PKGS$species_latin) %in% x)
	if(length(i)) {
		return(i[1])
	}

	# 9606
	i = which(as.character(BIOC_ANNO_PKGS$taxon_id) %in% x)
	if(length(i)) {
		return(i[1])
	}

	return(i)
}


# == title
# Extend TSS
#
# == param
# -orgdb Name of "org.*" packages from Bioconductor. All supported OrgDb packages are in ``rGREAT:::BIOC_ANNO_PKGS$orgdb``. 
# -verbose Whether to print messages.
# -... All pass to `extendTSS`.
#
# == value
# A `GenomicRanges::GRanges` object with one meta column 'gene_id'.
#
# == example
# if(FALSE) {
# extendTSSFromOrgDb("Org.Hs.eg.db")
# extendTSSFromOrgDb("hg19")
# }
extendTSSFromOrgDb = function(orgdb, verbose = great_opt$verbose, ...) {

	if(!inherits(orgdb, "character")) {
		orgdb = orgdb@.xData$packageName
	}

	i = detect_orgdb(orgdb)
	if(length(i) == 0) {
		stop_wrap(qq("OrgDb package '@{orgdb}' is not supported."))
	}
	orgdb = BIOC_ANNO_PKGS$orgdb[i]

	map_chr = get_table_from_orgdb("CHR$", orgdb)
	map_pos = get_table_from_orgdb("CHRLOC$", orgdb)
	chrlen = get_table_from_orgdb("CHRLENGTHS$", orgdb)

	if(orgdb == "org.Sc.sgd.db") {
		if(any(grepl("chrXIII", names(chrlen)))) {
			chr_nm_map = c("chrIV" = "chr4",
			               "chrXV" = "chr15",
			               "chrVII" = "chr7",
			               "chrXII" = "chr12",
			               "chrXVI" = "chr16",
			               "chrXIII" = "chr13",
			               "chrII" = "chr2", 
			               "chrXIV" = "chr14",
			               "chrX" = "chr10",
			               "chrXI" = "chr11",
			               "chrV" = "chr5",
			               "chrVIII" = "chr8",
			               "chrIX" = "chr9",
			               "chrIII" = "chr3",
			               "chrVI" = "chr6",
			               "chrI" = "chr1",
			               "chrM" = "chrM")
			names(chrlen) = chr_nm_map[ names(chrlen) ]
		}
	}

	if(is.null(map_chr) || is.null(map_pos) || is.null(chrlen)) {
		stop_wrap(qq("@{orgdb} does not contain gene coordinates information."))
	}

	hash = digest::digest(list(orgdb, ...))

	if(!is.null(rGREAT_env$extended_tss[[hash]])) {
		if(verbose) {
			message("* extended_tss is already cached, directly use it.")
		}
		return(rGREAT_env$extended_tss[[hash]])
	}

	if(!grepl("^chr", names(chrlen)[1])) {
		names(chrlen) = paste0("chr", names(chrlen))
	}
	
	map_pos = sapply(as.list(map_pos), function(x) {
		if(!is.null(names(x))) {
			x = x[nchar(names(x)) < 10]
		}
		x = x[x > 0]
		if(length(x)) {
			min(x)
		} else {
			NA
		}
	})
	map_pos = map_pos[!is.na(map_pos)]

	map_chr = as.list(map_chr)

	cn = intersect(names(map_pos), names(map_chr))
	map_pos = map_pos[cn]
	map_chr = sapply(map_chr[cn], function(x) x[1])

	if(!grepl("^chr", sample(map_chr, 1))) {
		map_chr = paste0("chr", map_chr)
	}
	
	strand = ifelse(map_pos > 0, "+", "-")
	map_pos = abs(map_pos)

	tss = GRanges(seqnames = map_chr, ranges = IRanges(map_pos, map_pos), strand = strand, gene_id = names(map_pos))

	l = map_pos <= chrlen[map_chr]; l[is.na(l)] = FALSE  # possible chr in map_chr but not in chrlen
	tss = tss[l]
	tss = unique(tss)

	chrlen = chrlen[nchar(names(chrlen)) < 10]

	chromosomes = names(chrlen)

	if(verbose) {
		message_str = qq("* restrict chromosomes to '@{paste(chromosomes, collapse=', ')}'.")
		message_str = strwrap(message_str, width = 90)
		if(length(message_str) > 1) {
			message_str[2:length(message_str)]= paste0("    ", message_str[2:length(message_str)])
		}
		message(paste(message_str, collapse = "\n"))
	}

	tss = tss[seqnames(tss) %in% chromosomes]

	map_genetype = get_table_from_orgdb("GENETYPE$", orgdb)
	if(!is.null(map_genetype)) {
		gene_type = unlist(as.list(map_genetype))
		ind = gene_type[tss$gene_id] == "protein-coding"
		ind[is.na(ind)] = FALSE
		tss = tss[ind]

		if(verbose) {
			message(qq("* @{length(tss)}/@{length(ind)} protein-coding genes left."))
		}
	}

	if(length(tss) == 0) {
		stop_wrap("No TSS left.")
	}

	if(verbose) {
		message(qq("* update seqinfo to the selected chromosomes."))
	}
	seqlevels(tss) = chromosomes
	info = seqinfo(tss)
	info = info[chromosomes, ]
	seqinfo(tss) = info

	seqlengths(tss) = chrlen

	extended_tss = extendTSS(tss, seqlengths(tss), verbose = verbose, 
		.attr = list(orgdb = orgdb, genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i], gene_id_type = "ENTREZ"), ...)

	rGREAT_env$extended_tss[[hash]] = extended_tss

	return(extended_tss)
}

# == title
# Extend TSS
#
# == param
# -df A bed-like data frame where the first three columns should be chromosomes, start positions, end positions.
#    It does not matter whether regions correspond to genes or TSS. 
# -seqlengths A named vector of chromosome lengths.
# -genome UCSC genome can be set here, then ``seqlengths`` will be automatically retrieved from UCSC server.
# -strand The strand information can be provided in ``df`` as a column named "strand" or as a column with "+"/"-"/"*", or the strand
#      information can be provided as a vector and be assigined to this argument.
# -gene_id The gene ID information can be provided in ``df`` as a column named "gene_id", or it can be provided as a vector and be assigned to this argument.
# -gene_id_type Gene ID types in ``df``. You need to set this argument if you use built-in gene sets in `great` so that genes can be correctly mapped.
#      The value can only be one of "SYMBOL", "ENTREZ", "ENSEMBL" and "REFSEQ".
# -verbose Whether to print messages.
# -... All pass to `extendTSS`.
#
# == value
# A `GenomicRanges::GRanges` object with one meta column 'gene_id'.
#
extendTSSFromDataFrame = function(df, seqlengths, genome = NULL, 
	strand = NULL, gene_id = NULL, 
	gene_id_type = NULL, verbose = great_opt$verbose, ...) {

	colnames(df) = tolower(colnames(df))
	
	df[, 1] = as.vector(df[, 1])
	if(!is.null(strand)) {
		df$strand = strand
	}
	if(!is.null(gene_id)) {
		df$gene_id = gene_id
	}

	if(missing(seqlengths)) {
		seqlengths = read.chromInfo(species = genome)$chr.len
		seqlengths = seqlengths[nchar(seqlengths) < 10]
	}

	df = df[ df[, 1] %in% names(seqlengths), , drop = FALSE]

	gene = GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df[, 3]))

	if(verbose) {
		message("* check strand column in df.")
	}
	if(!is.null(strand)) {
		strand(gene) = strand
	} else {
		if("strand" %in% colnames(df)) {
			strand(gene) = df[, "strand"]
		} else {
			l = !sapply(df, is.numeric)
			l[1] = FALSE
			flag = FALSE
			for(i in which(l)) {
				if(all(df[, i] %in% c("+", "-", "*"))) {
					flag = TRUE
					break
				}
			}
			if(flag) {
				strand(gene) = df[, i]
			} else {
				stop_wrap("Cannot find strand column in `df`.")
			}
		}
	}

	if(verbose) {
		message("* check gene_id column in df.")
	}
	
	if("gene_id" %in% colnames(df)) {
		gene$gene_id = df[, "gene_id"]
	} else {
		stop_wrap("Cannot find gene_id column in `df`")
	}
	
	gene$gene_id = as.character(gene$gene_id)
	names(gene) = gene$gene_id

	gene = gene[!is.na(names(gene))]

	if(verbose) {
		message("* add seqlengths information to the GRanges object.")
	}
	seqlevels(gene) = names(seqlengths)
	seqlengths(gene) = seqlengths

	if(!is.null(gene_id_type)) {
		if(!gene_id_type %in% c("SYMBOL", "ENTREZ", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`gene_id_type` can only be set as 'SYMBOL/ENTREZ/ENSEMBL/REFSEQ'")
		}
	}

	extendTSS(gene, seqlengths, verbose = verbose, .attr = list(genome = genome, gene_id_type = gene_id_type), ...)
}


# == title
# Extend TSS
#
# == param
# -gene A `GenomicRanges::GRanges` object of gene (or TSS) coordinates.
# -seqlengths A named vector of chromosome lengths. If it is not provided, it is taken by ``seqlengths(gene)``.
# -genome UCSC genome can be set here, then ``seqlengths`` will be automatically retrieved from UCSC server.
# -gene_id_type Gene ID types in ``gene``. You need to set this argument if you use built-in gene sets in `great` so that genes can be correctly mapped.
#      The value can only be one of "SYMBOL", "ENTREZ", "ENSEMBL" and "REFSEQ".
# -mode The mode to extend TSS. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See "Details" section.
# -basal_upstream In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.
# -basal_downstream In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.
# -extension Extensions from the basal domains. The value can also be a vector of length two which corresponds to extension to upstream and downstream respectively.
# -verbose Whether to print messages.
# -.attr Only used internally.
#
# == details
# Following are general explanations of the three modes for extending TSS:
#
# -``basalPlusExt`` 1. TSS are extended into basal domains (e.g. by upstream 5kb, downstream 1kb); 2. basal domains are sorted by their genomic
#      coordinates; 3. each basal domain is extended to its both sides until it reaches the next TSS's basal domain or it reaches the maximal extension (e.g. 1000kb).
# -``twoClosest`` 1. TSS are sorted by their genomic coordinates; 2. each TSS is extended to its both sides until it reaches the next TSS or it reaches the maximal extension (e.g. 1000kb).
# -``oneClosest`` 1. TSS are sorted by their genomic coordinates; 2. each TSS is extended to its both sides until it reaches the middle point of itself and the next TSS or it reaches the maximal extension (e.g. 1000kb).
#
# The official explanation is at https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules .
#
# == value
# A `GenomicRanges::GRanges` object with one meta column 'gene_id'.
#
extendTSS = function(gene, seqlengths = NULL, genome = NULL, 
	gene_id_type = NULL, mode = "basalPlusExt", basal_upstream = 5000, 
	basal_downstream = 1000, extension = 1000000,
	verbose = great_opt$verbose, .attr = list()) {

	if(is.null(seqlengths)) {
		if(is.null(genome)) {
			mt = attributes(metadata(gene))
			if(!is.null(mt$genome)) {
				genome = mt$genome
			}
		}

		if(is.null(genome)) {
			seqlengths = seqlengths(gene)
		} else {
			if(verbose) {
				message(qq("* use genome '@{genome}'."))
			}
			seqlengths = read.chromInfo(species = genome)$chr.len
			seqlengths = seqlengths[nchar(seqlengths) < 10]
		}

		# use the range of genes
		if(all(is.na(seqlengths))) {
			seqlengths = tapply(end(gene), seqnames(gene), max)
			seqlengths = structure(as.vector(seqlengths), names = names(seqlengths))
			message_wrap("! You havn't specified neither `seqlengths` nor `genome`, take the maximal ends of genes on chromosomes as the chromosome lengths.")
		}
	}
	sl = seqlengths

	if(is.null(gene_id_type)) {
		mt = attributes(metadata(gene))
		if(!is.null(mt$gene_id_type)) {
			gene_id_type = mt$gene_id_type
		}
	}

	if(!"gene_id" %in% colnames(mcols(gene))) {
		stop_wrap("`gene` must have a meta column `gene_id`")
	}
	names(gene) = gene$gene_id

	gene = gene[seqnames(gene) %in% names(sl)]

	tss = promoters(gene, upstream = 0, downstream = 1)

	if(verbose) {
		message(qq("* TSS extension mode is '@{mode}'."))
	}

	if(basal_downstream == 0) basal_downstream = 1
	if(mode == "basalPlusExt") {
		if(verbose) {
			message(qq("* construct the basal domains by extending @{basal_upstream}bp to upstream and @{basal_downstream}bp to downsteram of TSS."))
		}
		suppressWarnings(basal_tss <- promoters(tss, upstream = basal_upstream, downstream = basal_downstream))
		basal_tss = trim(basal_tss)
	} else if(mode == "twoClosest") {
		basal_tss = tss
	} else if(mode == "oneClosest") {
		basal_tss = tss
	} else {
		stop_wrap("`mode` should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'.")
	}

	strand(basal_tss) = "*"

	od = order(basal_tss)
	basal_tss = basal_tss[od]
	tss = tss[od]
	n = length(basal_tss)

	# distance to neighbours
	grl = split(basal_tss, seqnames(basal_tss))
	grl = grl[sapply(grl, length) > 0]

	if(verbose) {
		message("* calculate distances to neighbour regions.")
	}

	grl2 = lapply(grl, function(x) {
		names(x) = NULL
		as.data.frame(x)
	})
	gr_tb = NULL
	for(nm in names(grl2)) {
		gr = grl2[[nm]]
		n = nrow(gr)
		s = gr[, 2]
		e = gr[, 3]

		gr$dist_to_left = 0
		gr$dist_to_right = 0

		gr$dist_to_left[1] = s[1]
		gr$dist_to_right[n] = sl[nm] - e[n]
		if(n > 1) {
			if(mode == "oneClosest") {
				gr$dist_to_left[2:n] = floor((s[2:n] - e[1:(n-1)])/2)
				gr$dist_to_right[1:(n-1)] = floor((s[2:n] - e[1:(n-1)])/2)
			} else {
				gr$dist_to_left[2:n] = s[2:n] - e[1:(n-1)]
				gr$dist_to_right[1:(n-1)] = s[2:n] - e[1:(n-1)]
			}
		}

		gr$dist_to_left = ifelse(gr$dist_to_left > 0, gr$dist_to_left, 0)
		gr$dist_to_right = ifelse(gr$dist_to_right > 0, gr$dist_to_right, 0)

		gr_tb = rbind(gr_tb, gr)
	}

	basal_tss = as(gr_tb, "GRanges")

	if(length(extension) == 1) {
		extension = rep(extension, 2)
	}

	if(verbose) {
		message(qq("* extend to both sides until reaching the neighbour genes or to the maximal extension."))
	}
	extended_tss = basal_tss[, "gene_id"]

	extend_to_left = pmin(basal_tss$dist_to_left, ifelse(strand(tss) == "+", extension[1], extension[2]))
	start(extended_tss) = start(basal_tss) - extend_to_left

	extend_to_right = pmin(basal_tss$dist_to_right, ifelse(strand(tss) == "+", extension[2], extension[1]))
	end(extended_tss) = end(basal_tss) + extend_to_right

	names(extended_tss) = extended_tss$gene_id

	extended_tss$tss_position = start(tss)
	extended_tss$tss_strand = as.vector(strand(tss))

	seqlevels(extended_tss) = names(sl)
	suppressWarnings(seqlengths(extended_tss) <- sl)
	extended_tss = trim(extended_tss)

	if(!is.null(gene_id_type)) {
		if(!gene_id_type %in% c("SYMBOL", "ENTREZ", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`gene_id_type` can only be set as 'SYMBOL/ENTREZ/ENSEMBL/REFSEQ'")
		}
	}

	if(!is.null(gene_id_type)) {
		.attr$gene_id_type = gene_id_type
	}
	if(!is.null(genome)) {
		.attr$genome = genome
	}
	attributes(metadata(extended_tss)) = .attr

	attr(extended_tss, "mode") = list(
		"mode" = mode, 
		basal_upstream = as.integer(basal_upstream), 
		basal_downstream = as.integer(basal_downstream), 
		extension = as.integer(extension)
	)

	return(extended_tss)
}

