
# == title
# Get gap regions from UCSC
#
# == param
# -genome UCSC genome, such as "hg19".
# -seqnames A vector of chromosome names.
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == example
# getGapFromUCSC("hg19")
getGapFromUCSC = function(genome, seqnames = NULL) {

	lt = readRDS(system.file("extdata", "ucsc_gaps.rds", package = "rGREAT"))

	if(genome %in% names(lt)) {
		tb = lt[[genome]]
	} else {

		url = paste0("https://hgdownload.cse.ucsc.edu/goldenPath/", genome, "/database/gap.txt.gz")


		gap = paste0(tempdir(), "/", genome, "_gap.txt.gz")
		if(!file.exists(gap)) {
			e = try(suppressWarnings(download.file(url, destfile = gap, quiet = TRUE)), silent = TRUE)
			if(inherits(e, "try-error")) {
				stop_wrap(qq("It seems UCSC does not provide 'gap.txt.gz' for @{genome} or the internet connection was interrupted. If possible, download gap file directly from @{url}."))
			}
		}
		
		tb = read.table(gzfile(gap), sep = "\t", stringsAsFactors = FALSE)[, 2:4]
	}
	
	tb[, 2] = tb[, 2] + 1

	if(!is.null(seqnames)) {
		tb = tb[tb[, 1] %in% seqnames, , drop = FALSE]
	}
	GRanges(seqnames = tb[, 1], ranges = IRanges(tb[, 2], tb[, 3]))
}

getChromInfoFromUCSC = function(genome, seqlevels = NULL, max_seq = 500) {

	suppressWarnings(tb <- GenomeInfoDb::getChromInfoFromUCSC(genome, assembled.molecules.only = TRUE))
	if(genome %in% registered_UCSC_genomes()[, "genome"]) {
		chrlen = structure(tb[, 2], names = tb[, 1])
	} else {
		chrlen = structure(tb[, 2], names = tb[, 1])
		if(is.null(seqlevels)) {
			chrlen = filter_seqlength(chrlen, max_seq = max_seq)
		} else {
			chrlen = chrlen[intersect(names(chrlen), seqlevels)]
		}
	}
	chrlen
}

filter_seqlength = function(sl, max_seq = 500) {
	sl = sl*1.0
	sl = structure(sl, names = names(sl))

	sl = sort(sl)	
	l = cumsum(sl*1.0)/sum(sl*1.0) > 0.05 & nchar(names(sl)) < min(nchar(names(sl))) + 5
	if(sum(l) < 5) {
		l = cumsum(sl*1.0)/sum(sl*1.0) > 0.05
	}
	sl = sl[l]
	if(length(sl) > max_seq) {
		sl = sl[order(-sl)[1:max_seq]]
	}
	sl
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
getRefSeqGenesFromUCSC = function(genome, subset = c("RefSeqSelect", "RefSeqCurated")) {
	subset = match.arg(subset)[1]
	url = qq("https://hgdownload.cse.ucsc.edu/goldenPath/@{genome}/database/ncbi@{subset}.txt.gz")

	refseq = paste0(tempdir(), qq("/@{genome}_ncbi@{subset}.txt.gz"))
	if(!file.exists(refseq)) {
		e = try(suppressWarnings(download.file(url, destfile = refseq, quiet = TRUE)), silent = TRUE)
		if(inherits(e, "try-error")) {
			if(subset == "RefSeqSelect") {
				stop_wrap(qq("It seems UCSC does not provide 'ncbi@{subset}.txt.gz' for @{genome} or the internet connection was interrupted. If possible, download file directly from @{url}. Or you can try 'RefSeqCurated' subset."))
			} else {
				stop_wrap(qq("It seems UCSC does not provide 'ncbi@{subset}.txt.gz' for @{genome} or the internet connection was interrupted. If possible, download file directly from @{url}."))
			}
		}
	}
	
	tb = read.table(gzfile(refseq), sep = "\t", stringsAsFactors = FALSE)[, c(3, 5, 6, 4, 2)]
	tb[, 5] = gsub("\\.\\d+", "", tb[, 5])
	tb[, 2] = tb[, 2] + 1

	gene_id_type = "REFSEQ"

	i = which(BIOC_ANNO_PKGS$genome_version_in_txdb == genome)
	if(length(i)) {
		orgdb = BIOC_ANNO_PKGS$orgdb[i][1]
		map = get_table_from_orgdb("REFSEQ2EG$", orgdb)
		map = unlist(as.list(map))
		tb[, 5] = map[tb[, 5]]
		tb = tb[!is.na(tb[, 5]), ,drop = FALSE]
		gene_id_type = "ENTREZ"
	}
	
	gr = GRanges(seqnames = tb[, 1], ranges = IRanges(tb[, 2], tb[, 3]), strand = tb[, 4], gene_id = tb[, 5])
	gr = unique(gr)

	chrlen = getChromInfoFromUCSC(genome)

	gr = gr[seqnames(gr) %in% names(chrlen)]
	seqlevels(gr) = names(chrlen)
	seqlengths(gr) = chrlen

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

	tb = read.table(gzfile(system.file("extdata", paste0("GREATv4.genes.", genome, ".tsv.gz"), package = "rGREAT")))
	tb[is.na(tb[, 1]), 1] = ""
	if(genome %in% c("hg19", "hg38")) {
		map1 =  unlist(as.list(org.Hs.eg.db::org.Hs.egENSEMBL2EG))
		map2 =  unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG))
	} else {
		map1 =  unlist(as.list(org.Mm.eg.db::org.Mm.egENSEMBL2EG))
		map2 =  unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG))
	}

	gene_id_1 = map1[tb[, 1]]
	l = nchar(tb[, 1]) > 15 
	gene_id_1[l] = sapply(strsplit(tb[l, 1], ","), function(x) {
		x2 = map1[x]
		x2 = x2[!is.na(x2)]
		if(length(x2) == 0) {
			NA
		} else {
			x2[1]
		}
	})

	l = grep("^ENS.*\\d+$", tb[, 5])
	gene_id_2 = map2[tb[, 5]]
	gene_id_2[l] = map1[tb[l, 5]]

	gene_id = ifelse(is.na(gene_id_1), gene_id_2, gene_id_1)
	gene_id = unname(gene_id)

	gr = GRanges(seqnames = tb[, 2], ranges = IRanges(tb[, 3]+1, tb[, 3]+1), strand = tb[, 4], gene_id = gene_id)
	gr = gr[!is.na(gr$gene_id)]

	chrlen = getChromInfoFromUCSC(genome)

	gr = gr[seqnames(gr) %in% names(chrlen)]
	seqlevels(gr) = names(chrlen)
	seqlengths(gr) = chrlen

	m = metadata(gr)
	attr(m, "genome") = genome
	attr(m, "gene_id_type") = "ENTREZ"
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
getGenesFromGencode = function(version) {
	lt = readRDS(system.file("extdata", "gencode_gene.rds", package = "rGREAT"))

	if(version == "hg19") {
		versio = "v19"
	} else if(version == "hg38") {
		version = "v40"
	} else if(version == "mm9") {
		version = "vM1"
	} else if(version == "mm10") {
		version = "vM25"
	} else if(version == "mm11") {
		version = "vM29"
	}

	version = gsub("^V", "v", version)

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

	chrlen = getChromInfoFromUCSC(genome)

	gr = gr[seqnames(gr) %in% names(chrlen)]
	seqlevels(gr) = names(chrlen)
	seqlengths(gr) = chrlen

	m = metadata(gr)
	attr(m, "genome") = genome
	attr(m, "gene_id_type") = "ENSEMBL"
	metadata(gr) = m
	gr
}


getTSSFromOrgDb = function(orgdb) {
	stop_wrap("OrgDb is not supported any more.")
	
	i = detect_orgdb(orgdb)
	if(length(i) == 0) {
		stop_wrap(qq("OrgDb package '@{orgdb}' is not supported."))
	}
	orgdb = BIOC_ANNO_PKGS$orgdb[i]

	map_chr = get_table_from_orgdb("CHR$", orgdb)
	map_pos_start = get_table_from_orgdb("CHRLOC$", orgdb)
	map_pos_end = get_table_from_orgdb("CHRLOCEND$", orgdb)
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

	if(is.null(map_chr) || is.null(map_pos_start) || is.null(chrlen)) {
		stop_wrap(qq("@{orgdb} does not contain gene coordinates information."))
	}

	if(!grepl("^chr", names(chrlen)[1])) {
		names(chrlen) = paste0("chr", names(chrlen))
	}
	
	map_pos_start = sapply(as.list(map_pos_start), function(x) {
		if(!is.null(names(x))) {
			x = x[nchar(names(x)) < 10]
		}
		if(length(x)) {
			min(x)
		} else {
			NA
		}
	})
	map_pos_end = sapply(as.list(map_pos_end), function(x) {
		if(!is.null(names(x))) {
			x = x[nchar(names(x)) < 10]
		}
		if(length(x)) {
			max(x)
		} else {
			NA
		}
	})
	l = !is.na(map_pos_start) & !is.na(map_pos_end)
	map_pos_start = map_pos_start[l]
	map_pos_end = map_pos_end[l]

	map_chr = as.list(map_chr)

	cn = intersect(names(map_pos_start), names(map_chr))
	map_pos_start = map_pos_start[cn]
	map_pos_end = map_pos_end[cn]
	map_chr = sapply(map_chr[cn], function(x) x[1])

	if(!grepl("^chr", sample(map_chr, 1))) {
		map_chr = paste0("chr", map_chr)
	}
	
	strand = ifelse(map_pos_start > 0, "+", "-")

	gene = GRanges(seqnames = map_chr, ranges = IRanges(abs(map_pos_start), abs(map_pos_end)), strand = strand, gene_id = names(map_pos_start))
	tss = promoters(gene, upstream = 0, downstream = 1)

	l = map_pos_end <= chrlen[map_chr]; l[is.na(l)] = FALSE  # possible chr in map_chr but not in chrlen
	tss = tss[l]
	tss = unique(tss)

	chrlen = chrlen[nchar(names(chrlen)) < 10]

	chromosomes = names(chrlen)

	tss = tss[seqnames(tss) %in% chromosomes]

	map_genetype = get_table_from_orgdb("GENETYPE$", orgdb)
	if(!is.null(map_genetype)) {
		gene_type = unlist(as.list(map_genetype))
		ind = gene_type[tss$gene_id] == "protein-coding"
		ind[is.na(ind)] = FALSE
		tss = tss[ind]
	}

	if(length(tss) == 0) {
		stop_wrap("No TSS left.")
	}

	seqlevels(tss) = chromosomes
	info = seqinfo(tss)
	info = info[chromosomes, ]
	seqinfo(tss) = info

	seqlengths(tss) = chrlen

	tss
}

getTSSFromTxDb = function(txdb_pkg) {

	i = detect_txdb(txdb_pkg)
	if(length(i) == 0) {
		stop_wrap(qq("TxDb package '@{txdb_pkg}' is not supported."))
	}

	txdb_pkg = BIOC_ANNO_PKGS$txdb[i]

	check_pkg(txdb_pkg, bioc = TRUE)

	suppressMessages(gene <- genes( getFromNamespace(txdb_pkg, txdb_pkg) ))

	genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i]
	cl = CHR_LEN_DB[[genome]]
	chromosomes = names(cl)

	gene = gene[seqnames(gene) %in% chromosomes]

	orgdb = BIOC_ANNO_PKGS$orgdb[i]

	check_pkg(orgdb, bioc = TRUE)

	all_tb = ls(envir = getNamespace(orgdb))
	i_tb = which(grepl("ENSEMBL2EG", all_tb))

	gene_id_type = BIOC_ANNO_PKGS$gene_id_in_txdb[i]
	if(BIOC_ANNO_PKGS$gene_id_in_txdb[i] == "Ensembl gene ID" & length(i_tb)) {


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

	}


	seqlevels(gene) = chromosomes
	info = seqinfo(gene)
	info = info[chromosomes, ]
	seqinfo(gene) = info

	promoters(gene, upstream = 0, downstream = 1)
}


getGenesFromBioMart = function(dataset, filter = FALSE, max_seq = 500) {
	check_pkg("BioMartGOGeneSets", bioc = TRUE, github = "jokergoo")
	g = BioMartGOGeneSets::getBioMartGenes(dataset)
	colnames(mcols(g))[1] = "gene_id"
	if(filter) {
		sl = tapply(end(g), seqnames(g), max)
		sl = filter_seqlength(sl, max_seq = max_seq)
		g = g[seqnames(g) %in% names(sl)]

		seqlevels(g) = unique(as.character(seqnames(g)))
	}
	g
}


# == title
# Get gene sets from BioMart
#
# == param
# -dataset Name of the dataset.
# -ontology Value should be bp, mf or cc.
#
# == details
# GO gene sets are from ``BioMartGOGeneSets::getBioMartGOGeneSets``.
#
# == value
# A list of vectors where each vector contains Ensembl IDs annotated to a GO term.
getGeneSetsFromBioMart = function(dataset, ontology = "bp") {
	check_pkg("BioMartGOGeneSets", bioc = TRUE, github = "jokergoo")

	ontology = tolower(ontology)
	BioMartGOGeneSets::getBioMartGOGeneSets(dataset, ontology)
}

# == title
# Get the internally used TSS
#
# == param
# -tss_source The same format as in `great`.
# -biomart_dataset The same format as in `great`.
#
# == value
# A `GenomicRanges::GRanges` object.
#
getTSS = function(tss_source, biomart_dataset = NULL) {
	if(!is.null(biomart_dataset)) {
		
		biomart_dataset = tolower(biomart_dataset)
		genes = getGenesFromBioMart(biomart_dataset)
		tss = promoters(genes, upstream = 0, downstream = 1)
	} else {

		tss_source = parse_tss_source(tss_source)

		if(tss_source$category == "TxDb") {
			tss = getTSSFromTxDb(tss_source$source)

		} else if(tss_source$category == "OrgDb") {

			tss = getTSSFromOrgDb(tss_source$source)

		} else if(tss_source$category == "Gencode") {

			genes = getGenesFromGencode(tss_source$source)
			tss = promoters(genes, upstream = 0, downstream = 1)

		} else if(tss_source$category == "RefSeq") {

			genes = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqSelect")
			tss = promoters(genes, upstream = 0, downstream = 1)
			
		} else if(tss_source$category == "RefSeqCurated") {

			genes = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqCurated")
			tss = promoters(genes, upstream = 0, downstream = 1)

		} else if(tss_source$category == "RefSeqSelect") {

			genes = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqSelect")
			tss = promoters(genes, upstream = 0, downstream = 1)

		} else if(tss_source$category == "GREAT") {

			tss = getGREATDefaultTSS(tss_source$genome)
			
		} else {
			stop_wrap("Wrong `tss_source`.")
		}
	}
	tss
}
