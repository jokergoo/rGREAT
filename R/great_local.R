

if(identical(topenv(), .GlobalEnv)) {
	BIOC_ANNO_PKGS = read.table("~/project/development/rGREAT/inst/extdata/bioc_anno_pkgs.csv", header = TRUE, sep = "\t")
} else {
	BIOC_ANNO_PKGS = read.table(system.file("extdata", "bioc_anno_pkgs.csv", package = "rGREAT"), header = TRUE, sep = "\t")
}

rGREAT_env$extended_tss = list()


detect_txdb = function(x) {
	i = which(BIOC_ANNO_PKGS$txdb %in% x)
	if(length(i)) {
		return(i)
	}

	i = which(BIOC_ANNO_PKGS$genome_version_in_txdb %in% x)
	if(length(i)) {
		return(i[1])
	}

	i = which(grepl(x, BIOC_ANNO_PKGS$txdb))
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
# -chromosomes The original TxDb data contains not only the "normal chromosomes" but also a lot of small contigs. This
#         function can automatically identify those "normal chromosomes" (not guaranteed), but users can always control which chromosomes to use by specify a vector of chromosome names to this argument.
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
extendTSSFromTxDb = function(txdb, chromosomes = NULL, 
	verbose = great_opt$verbose, ...) {

	txdb_pkg = txdb
	if(!inherits(txdb_pkg, "character")) {
		txdb_pkg = txdb_pkg$packageName
	}

	i = detect_txdb(txdb_pkg)
	if(length(i) == 0) {
		stop_wrap(qq("TxDb package '@{txdb_pkg}' is not supported."))
	}

	txdb_pkg = BIOC_ANNO_PKGS$txdb[i]

	hash = digest::digest(list(txdb_pkg, chromosomes, ...))

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
	sl = seqlengths(gene)

	if(is.null(chromosomes)) {
		chromosomes = intersect(names(sl), guess_major_chromosomes(sl))
	} else {
		chromosomes = intersect(chromosomes, names(sl))
	}
	if(verbose) {
		message_str = qq("* restrict chromosomes to '@{paste(chromosomes, collapse=', ')}'.")
		message_str = strwrap(message_str)
		if(length(message_str) > 1) {
			message_str[2:length(message_str)]= paste0("    ", message_str[2:length(message_str)])
		}
		message(paste(message_str, collapse = "\n"))
	}

	gene = gene[seqnames(gene) %in% chromosomes]

	org_db = BIOC_ANNO_PKGS$orgdb[i]

	check_pkg(org_db, bioc = TRUE)

	if(BIOC_ANNO_PKGS$gene_id_in_txdb[i] == "Ensembl gene ID") {
		if(verbose) {
			message("* convert gene ID from 'Ensembl gene ID' to 'Entrez gene ID'.")
		}
		map = get_table_from_org_db("ENSEMBL2EG", org_db)
		map = sapply(map, function(x) x[1])

		gene_id = names(gene)
		new_gene_id = map[gene_id]

		l = !is.na(new_gene_id)

		gene = gene[l, ]
		gene$gene_id = new_gene_id[l]
		names(gene) = new_gene_id[l]
	}

	if(verbose) {
		message(qq("* check whether '*GENETYPE' table is available in '@{org_db}'."))
	}
	all_tb = ls(envir = getNamespace(org_db))
	i_tb = which(grepl("GENETYPE", all_tb))
	if(length(i_tb)) {
		gene_type_tb = getFromNamespace(all_tb[i_tb], org_db)
		gene_type = unlist(as.list(gene_type_tb))
		ind = gene_type[names(gene)] == "protein-coding"
		ind[is.na(ind)] = FALSE
		gene = gene[ind]

		if(verbose) {
			message(qq("* @{length(gene)} protein-coding genes left."))
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

	metadata(gene) = BIOC_ANNO_PKGS[i, , drop = FALSE]

	extended_tss = extendTSS(gene, seqlengths(gene), verbose = verbose, ...)

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
# -strand The strand information can be provided in ``df`` as a column named "strand" or as a column with "+"/"-"/"*", or the strand
#      information can be provided as a vector and be assigined to this argument.
# -gene_id The gene ID information can be provided in ``df`` as a column named "gene_id", or it can be provided as a vector and be assigned to this argument.
# -verbose Whether to print messages.
# -... All pass to `extendTSS`.
#
# == value
# A `GenomicRanges::GRanges` object with one meta column 'gene_id'.
#
extendTSSFromDataFrame = function(df, seqlengths, strand = NULL, gene_id = NULL, 
	verbose = great_opt$verbose, ...) {

	df[, 1] = as.vector(df[, 1])

	if(length(setdiff(df[, 1], names(seqlengths)))) {
		stop_wrap("All chromsomes in `df` should be included in `seqlengths`.")
	}
	gene = GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df[, 3]))

	colnames(df) = tolower(colnames(df))

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
	if(!is.null(gene_id)) {
		gene$gene_id = gene_id
	} else {
		if("gene_id" %in% colnames(df)) {
			gene$gene_id = df[, "gene_id"]
		} else {
			stop_wrap("Cannot find gene_id column in `df`")
		}
	}

	gene$gene_id = as.character(gene$gene_id)
	names(gene) = gene$gene_id

	if(verbose) {
		message("* add seqlengths infomation to the GRanges object.")
	}
	seqlevels(gene) = names(seqlengths)
	seqlengths(gene) = seqlengths

	extendTSS(gene, seqlengths, verbose = verbose, ...)
}


# == title
# Extend TSS
#
# == param
# -gene A `GenomicRanges::GRanges` object of gene (or TSS) coordinates.
# -seqlengths A named vector of chromosome lengths. If it is not provided, it is taken by ``seqlengths(gene)``.
# -mode The mode to extend TSS. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See "Details" section.
# -basal_upstream In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.
# -basal_downstream In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.
# -extension Extensions from the basal domains.
# -verbose Whether to print messages.
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
extendTSS = function(gene, seqlengths, mode = "basalPlusExt",
	basal_upstream = 5000, basal_downstream = 1000, extension = 1000000,
	verbose = great_opt$verbose) {

	if(missing(seqlengths)) {
		seqlengths = seqlengths(gene)
	}
	sl = seqlengths

	if(!"gene_id" %in% colnames(mcols(gene))) {
		stop_wrap("`gene` must have a meta column `gene_id`")
	}
	names(gene) = gene$gene_id

	tss = promoters(gene, upstream = 0, downstream = 1)

	if(verbose) {
		message(qq("* TSS extension mode is '@{mode}'."))
	}

	if(mode == "basalPlusExt") {
		if(verbose) {
			message(qq("* construct the basal domains by extending @{basal_upstream}bp to upstream and @{basal_downstream}bp to downsteram of TSS."))
		}
		basal_tss = promoters(tss, upstream = basal_upstream, downstream = basal_downstream)
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

	if(length(setdiff(names(grl), names(sl)))) {
		stop_wrap("`seqlengths` should cover all chromosomes in `gene`.")
	}

	if(verbose) {
		message("* calculate distances to neighbour regions.")
	}
	for(nm in names(grl)) {
		gr = grl[[nm]]
		n = length(gr)
		s = start(gr)
		e = end(gr)

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
		grl[[nm]] = gr
	}

	basal_tss = unname(unlist(grl))

	if(verbose) {
		message(qq("* extend to both sides until reaching the neighbour genes or to the maximal extension (@{as.integer(extension)}bp)."))
	}
	extended_tss = basal_tss[, "gene_id"]

	extend_to_left = pmin(basal_tss$dist_to_left, extension)
	start(extended_tss) = start(basal_tss) - extend_to_left + 1

	extend_to_right = pmin(basal_tss$dist_to_right, extension)
	end(extended_tss) = end(basal_tss) + extend_to_right - 1

	names(extended_tss) = extended_tss$gene_id

	extended_tss$tss_position = start(tss)
	extended_tss$tss_strand = as.vector(strand(tss))

	seqlevels(extended_tss) = names(sl)
	seqlengths(extended_tss) = sl

	attr(extended_tss, "mode") = list(
		"mode" = mode, 
		basal_upstream = as.integer(basal_upstream), 
		basal_downstream = as.integer(basal_downstream), 
		extension = as.integer(extension)
	)

	return(extended_tss)
}

# sl: a vector of chromosome length (the vector with names)
guess_major_chromosomes = function(sl) {
	sl = sort(sl, decreasing = TRUE)
	n = length(sl)
	if(n == 1) {
		return(names(sl))
	}

	ind = which(sl[1:(n-1)]/sl[2:n] > 5)
	names(sl)[1:ind[1]]

}


# == title
# Perform GREAT analysis
#
# == param
# -gr A `GenomicRanges::GRanges` object. This is the input regions.
# -gene_sets A single string of defautly supported gene sets collections (see the full list in "Details" section), or a named list of vectors where each vector correspond to a gene set.
# -txdb Name of "TxDb.*" packages from Bioconductor. All supported TxDb packages are in ``rGREAT:::BIOC_ANNO_PKGS$txdb``. Note short genome version can also be used
#     here such as "hg19" or "hg19.knownGene".
# -chromosomes The original TxDb data contains not only the "normal chromosomes" but also a lot of small contigs. This
#         function can automatically identify those "normal chromosomes" (not guaranteed), but users can always control which chromosomes to use by specify a vector of chromosome names to this argument.
# -mode The mode to extend genes. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See `extendTSS` for details.
# -basal_upstream In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.
# -basal_downstream In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.
# -extension Extensions from the basal domains.
# -extended_tss If your organism is not defaultly supported, you can first prepare one by `extendTSSFromDataFrame` or `extendTSS`,
#      and set the object to this argument.
# -background Background regions.
# -exclude Regions that are excluded from analysis such as gap regions (which can be get by `getGapFromUCSC`).
# -verbose Whether to print messages.
#
# == details
# When ``background`` or ``exclude`` is set, the analysis is restricted in the background regions, still by using Binomial method. Note
# this is different from the original GREAT method which uses Fisher's exact test if background regions is set. See `submitGreatJob` for explanations.
#
# rGREAT supports the following GO gene sets for all organisms (note "GO:" can be omitted):
#
# -"GO:BP": Biological Process, from GO.db package.
# -"GO:CC": Cellular Component, from GO.db package.
# -"GO:MP": Molecular Function, from GO.db pacakge.
#
# rGREAT also supports gene sets collections from MSigDB (with the msigdbr package, note this is only for human, "msigdb:" can be omitted):
#
# -"msigdb:H" Hallmark gene sets.
# -"msigdb:C1" Positional gene sets.
# -"msigdb:C2" Curated gene sets.
# -"msigdb:C2:CGP" C2 subcategory: chemical and genetic perturbations gene sets. 
# -"msigdb:C2:CP" C2 subcategory: canonical pathways gene sets. 
# -"msigdb:C2:CP:BIOCARTA" C2 subcategory: BioCarta subset of CP.
# -"msigdb:C2:CP:KEGG" C2 subcategory: KEGG subset of CP.
# -"msigdb:C2:CP:PID" C2 subcategory: PID subset of CP.
# -"msigdb:C2:CP:REACTOME" C2 subcategory: REACTOME subset of CP.
# -"msigdb:C2:CP:WIKIPATHWAYS" C2 subcategory: WIKIPATHWAYS subset of CP.
# -"msigdb:C3" Regulatory target gene sets.
# -"msigdb:C3:MIR:MIRDB" miRDB of microRNA targets gene sets.
# -"msigdb:C3:MIR:MIR_LEGACY" MIR_Legacy of MIRDB.
# -"msigdb:C3:TFT:GTRD" GTRD transcription factor targets gene sets.
# -"msigdb:C3:TFT:TFT_LEGACY" TFT_Legacy.
# -"msigdb:C4" Computational gene sets.
# -"msigdb:C4:CGN" C4 subcategory: cancer gene neighborhoods gene sets.
# -"msigdb:C4:CM" C4 subcategory: cancer modules gene sets.
# -"msigdb:C5" Ontology gene sets.
# -"msigdb:C5:GO:BP" C5 subcategory: BP subset.
# -"msigdb:C5:GO:CC" C5 subcategory: CC subset.
# -"msigdb:C5:GO:MF" C5 subcategory: MF subset.
# -"msigdb:C5:HPO" C5 subcategory: human phenotype ontology gene sets.
# -"msigdb:C6" Oncogenic signature gene sets.
# -"msigdb:C7" Immunologic signature gene sets.
# -"msigdb:C7:IMMUNESIGDB" ImmuneSigDB subset of C7.
# -"msigdb:C7:VAX" C7 subcategory: vaccine response gene sets.
# -"msigdb:C8" Cell type signature gene sets.
#
# If the defaultly supported gene sets and TxDb are used, Entrez gene ID is always used as the main gene ID. If you provide a self-defined
# ``gene_sets`` or ``extended_tss``, you need to make sure they two have the same gene ID types.
#
# == value
# A `GreatObject-class` object. The following methods can be applied on it:
# 
# - `getEnrichmentTable,GreatObject-method` to retrieve the result table. 
# - `getRegionGeneAssociations,GreatObject-method` to get the associations between input regions and genes.
# - `plotRegionGeneAssociations,GreatObject-method` to plot the associations bewteen input regions and genes.
# - `shinyReport,GreatObject-method` to view the results by a shiny application.
#
# == example
# if(FALSE) {
# gr = randomRegions()
# res = great(gr, "MSigDB:H", "hg19")
# }
great = function(gr, gene_sets, txdb, chromosomes = NULL,
	mode = "basalPlusExt", basal_upstream = 5000, 
	basal_downstream = 1000, extension = 1000000,
	extended_tss = NULL, background = NULL, exclude = NULL,
	verbose = great_opt$verbose) {

	if(is.null(extended_tss)) {
		extended_tss = extendTSSFromTxDb(txdb, chromosomes = chromosomes, mode = mode, basal_upstream = basal_upstream,
			basal_downstream = basal_downstream, extension = extension)
	}

	if(is.character(gene_sets) && is.atomic(gene_sets)) {

		gene_sets_name = gene_sets
		mtb = metadata(extended_tss)

		gene_sets = get_defaultly_suppported_gene_sets(gene_sets, mtb, verbose = verbose)

		if(is.null(gene_sets)) {
			stop_wrap("Please set `gene_sets` as a named list of vectors where each vector corresponds to a gene set. The gene ID type should be the same as in `txdb` or `extended_tss`. If the defaultly supported TxDb package is used, Entrez gene ID is the main gene ID.")
		}
		
	} else {
		gene_sets_name = "self-provided"

		if("gene_id" %in% colnames(mcols(extended_tss))) {
			stop_wrap("`extended_tss` must have a meta column `gene_id`")
		}
		names(extended_tss) = extended_tss$gene_id

		if(verbose) {
			message("* check gene ID type in `gene_sets` and in `extended_tss`.")
		}

		v1 = unique(unlist(gene_sets[seq_len(min(10, length(gene_sets)))]))
		if(length(intersect(v1, extended_tss$gene_id))/length(v1) < 0.2) {
			warning_wrap("It seems the gene ID type in `gene_sets` is different from in `extended_tss`.")
		}
	}

	mode_param = attr(extended_tss, "mode")

	sl = seqlengths(extended_tss)
	if(any(is.na(sl))) {
		stop_wrap("`extended_tss` must have `seqlengths` been set.")
	}

	gr_genome = GRanges(seqnames = names(sl), ranges = IRanges(1, sl))
	if(is.null(background)) {
		background = gr_genome

		if(verbose) {
			message("* use whole genome as background.")
		}
	} else {
		if(verbose) {
			message("* use self-defined background regions.")
		}
	}
	strand(background) = "*"
	if(!is.null(exclude)) {
		strand(exclude) = "*"
		background = setdiff(background, exclude)

		if(verbose){
			message("* remove excluded regions from background.")
		}
	}

	background_total_length = sum(width(background))

	strand(gr) = "*"

	gr_origin = gr

	gr_mid = gr
	start(gr_mid) = end(gr_mid) = mid(gr)
	gr = gr_mid

	if(verbose) {
		message("* overlap `gr` to background regions.")
	}
	gr = intersect(gr, background)

	n_total = length(gr)
	
	if(verbose) {
		message("* check which genes are in the gene sets.")
	}
	gene_sets = lapply(gene_sets, as.character)
	gene_sets_ind = lapply(gene_sets, function(x) which(names(extended_tss) %in% x))

	l = sapply(gene_sets_ind, function(x) length(x) > 0)
	gene_sets_ind = gene_sets_ind[l]
	gene_sets = lapply(gene_sets_ind, function(x) names(extended_tss)[x])

	n_set = length(gene_sets)
	
	if(verbose) {
		message(qq("* in total @{n_set} gene sets."))
	}

	p = numeric(n_set)
	n_obs = numeric(n_set)
	n_exp = numeric(n_set)
	fold = numeric(n_set)
	prop = numeric(n_set)

	gene_hits = numeric(n_set)


	if(verbose) {
		message("* overlap extended TSS to background regions.")
	}
	ov2 = findOverlaps(extended_tss, background)
	r1 = extended_tss[queryHits(ov2)]
	r2 = background[subjectHits(ov2)]

	extended_tss2 = pintersect(r1, r2)
	extended_tss2$gene_id = r1$gene_id
	names(extended_tss2) = r1$gene_id
	extended_tss = extended_tss2

	if(verbose) {
		message(qq("* overlap `gr` to every extended TSS."))
	}
	ov = findOverlaps(gr, extended_tss)

	if(verbose) {
		message("* perform Binomial test for each biological term.")
	}
	pb = progress::progress_bar$new(total = n_set)
	for(i in 1:n_set) {

		if(verbose) pb$tick()

		ind = gene_sets_ind[[i]]
		fgr = extended_tss[ind]
		fgr = reduce(fgr)
		prop[i] = sum(width(fgr))/background_total_length

		n_hits = length(unique(queryHits(ov)[subjectHits(ov) %in% ind]))

		if(n_hits == 0) {
			p[i] = 1
		} else {
			p[i] = 1 - pbinom(n_hits - 1, n_total, prop[i])
		}
		n_obs[i] = n_hits
		n_exp[i] = prop[i]*n_total

		gene_hits[i] = length(intersect(subjectHits(ov), ind))
	}

	fold = n_obs/n_exp

	df = data.frame(id = names(gene_sets),
			   genome_fraction = prop,
		       observed_region_hits = n_obs,
		       fold_enrichment = fold,
		       p_value = p,
		       p_adjust = p.adjust(p, "BH"),
		       observed_gene_hits = gene_hits,
		       gene_set_size = sapply(gene_sets_ind, length))

	df = df[order(df$p_adjust, df$p_value, -df$fold_enrichment), , drop = FALSE]

	if(all(grepl("^GO:\\d+$", df$id))) {
		check_pkg("GO.db", bioc = TRUE)
		terms = Term(GO.db::GOTERM)

		df = cbind(data.frame(id = df$id, description = terms[df$id]),
			       df[, -1])
	}
	rownames(df) = NULL

	mtb = metadata(extended_tss)
	if(length(mtb) == 0) {
		txdb = "unknown"
		orgdb = "unknown"
	} else if(!is.data.frame(mtb)) {
		txdb = "unknown"
		orgdb = "unknown"
	} else if(!"orgdb" %in% names(mtb)) {
		txdb = "unknown"
		orgdb = "unknown"
	} else {
		txdb = mtb$txdb
		orgdb = mtb$orgdb
	}

	obj = GreatObject(
		table = df,
		gr = gr_origin,
		gene_sets = gene_sets,
		gene_sets_name = gene_sets_name,
		extended_tss = extended_tss,  # has been intersected with background
		background = background,
		txdb = txdb,
		orgdb = orgdb,
		param = mode_param
	)

	return(obj)
}

# support two types
# -GO (GO:BP, BP, GO:CC, CC, GO:MF, MF)
# -MsigDB
get_defaultly_suppported_gene_sets = function(name, mtb, verbose = great_opt$verbose) {
	if(!is.data.frame(mtb)) {
		return(NULL)
	}
	if(!"orgdb" %in% names(mtb)) {
		return(NULL)
	}
	name = tolower(name)

	if(grepl("^c\\d", name) || name == "h") {
		name = paste0("msigdb:", name)
	}

	if(name %in% c("go:bp", "bp", "go:cc", "cc", "go:mf", "mf")) {
		lt = as.list(get_table_from_org_db("GO2ALL.*S$", mtb$orgdb))

		check_pkg("GO.db", bioc = TRUE)

		ontology = Ontology(GO.db::GOTERM)

		if(name %in% c("go:bp", "bp")) {
			l = ontology[names(lt)] == "BP"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:BP ontology with @{length(gene_sets)} gene sets (source: @{mtb$orgdb})."))
			}
		} else if(name %in% c("go:mf", "mf")) {
			l = ontology[names(lt)] == "MF"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:MF ontology with @{length(gene_sets)} gene sets (source: @{mtb$orgdb})."))
			}
		} else if(name %in% c("go:cc", "cc")) {
			l = ontology[names(lt)] == "CC"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:CC ontology with @{length(gene_sets)} gene sets (source: @{mtb$orgdb})."))
			}
		}

		return(gene_sets)
	} else if(grepl("msigdb", name)) {
		if(mtb$species_name != "human") {
			stop_wrap("MSigDB only supports human.")
		}

		check_pkg("msigdbr", bioc = FALSE)

		if(name == "msigdb:h") {
			tb = msigdbr(category = "H")
		} else if(name == "msigdb:c2") {
			tb = msigdbr(category = "C2")
		} else if(name == "msigdb:c2:cgp") {
			tb = msigdbr(category = "C2", subcategory = "CGP")
		} else if(name == "msigdb:c2:cp") {
			tb = msigdbr(category = "C2", subcategory = "CP")
		} else if(name == "msigdb:c2:cp:biocarta" || name == "msigdb:c2:biocarta") {
			tb = msigdbr(category = "C2", subcategory = "CP:BIOCARTA")
		} else if(name == "msigdb:c2:cp:kegg" || name == "msigdb:c2:cp:kegg") {
			tb = msigdbr(category = "C2", subcategory = "CP:KEGG")
		} else if(name == "msigdb:c2:cp:pid" || name == "msigdb:c2:cp:pid") {
			tb = msigdbr(category = "C2", subcategory = "CP:PID")
		} else if(name == "msigdb:c2:cp:reactome" || name == "msigdb:c2:cp:reactome") {
			tb = msigdbr(category = "C2", subcategory = "CP:REACTOME")
		} else if(name == "msigdb:c2:cp:wikipathways" || name == "msigdb:c2:cp:wikipathways") {
			tb = msigdbr(category = "C2", subcategory = "CP:WIKIPATHWAYS")
		} else if(name == "msigdb:c3") {
			tb = msigdbr(category = "C3")
		} else if(name == "msigdb:c3:mir:mirdb") {
			tb = msigdbr(category = "C3", subcategory = "MIR:MIRDB")
		} else if(name == "msigdb:c3:mir:mir_legacy") {
			tb = msigdbr(category = "C3", subcategory = "MIR:MIR_Legacy")
		} else if(name == "msigdb:c3:tft:gtrd") {
			tb = msigdbr(category = "C3", subcategory = "TFT:GTRD")
		} else if(name == "msigdb:c3:tft:tft_legacy") {
			tb = msigdbr(category = "C3", subcategory = "TFT:TFT_Legacy")
		} else if(name == "msigdb:c4") {
			tb = msigdbr(category = "C4")
		} else if(name == "msigdb:c4:cgn") {
			tb = msigdbr(category = "C4", subcategory = "CGN")
		} else if(name == "msigdb:c4:cm") {
			tb = msigdbr(category = "C4", subcategory = "CM")
		} else if(name == "msigdb:c5") {
			tb = msigdbr(category = "C5")
		} else if(name == "msigdb:c5:go:bp") {
			tb = msigdbr(category = "C5", subcategory = "GO:BP")
		} else if(name == "msigdb:c5:go:cc") {
			tb = msigdbr(category = "C5", subcategory = "GO:CC")
		} else if(name == "msigdb:c5:go:mf") {
			tb = msigdbr(category = "C5", subcategory = "GO:MF")
		} else if(name == "msigdb:c5:hpo") {
			tb = msigdbr(category = "C5", subcategory = "hpo")
		} else if(name == "msigdb:c6") {
			tb = msigdbr(category = "C6")
		} else if(name == "msigdb:c7") {
			tb = msigdbr(category = "C7")
		} else if(name == "msigdb:c7:immunesigdb") {
			tb = msigdbr(category = "C7", subcategory = "IMMUNESIGDB")
		} else if(name == "msigdb:c7:vax") {
			tb = msigdbr(category = "C7", subcategory = "VAX")
		} else if(name == "msigdb:c8") {
			tb = msigdbr(category = "C8")
		} else {
			stop_wrap("Wrong MSigDB gene set collection.")
		}

		gene_sets = split(tb$entrez_gene, tb$gs_name)

		if(verbose) {
			message(qq("* use @{name} collection with @{length(gene_sets)} gene sets (source: msigdbr)."))
		}
		return(gene_sets)
	}
}

# == title
# Class for local GREAT analysis
#
# == details
# `great` returns A `GreatObject-class` object.
#
GreatObject = setClass("GreatObject",
    slots = list(table = "data.frame",
    	         gr = "GRanges",
    	         gene_sets = "list",
    	         extended_tss = "GRanges",
    	         gene_sets_name = "character",
    	         background = "GRanges",
    	         txdb = "character",
    	         orgdb = "character",
    	         param = "ANY")
)


# == title
# Constructor method for GreatObject class
#
# == param
# -... arguments.
#
# == details
# There are following methods that can be applied on `GreatObject-class` object:
#
# - `getEnrichmentTable,GreatObject-method` to retrieve the result table. 
# - `getRegionGeneAssociations,GreatObject-method` to get the associations between input regions and genes.
# - `plotRegionGeneAssociations,GreatObject-method` to plot the associations bewteen input regions and genes.
# - `shinyReport,GreatObject-method` to view the results by a shiny application.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
GreatObject = function(...) {
    new("GreatObject", ...)
}

setMethod(f = "show",
	signature = "GreatObject",
	definition = function(object) {

	qqcat("Associate @{length(object@gr)} regions to @{length(object@gene_sets)} gene sets.\n")
	qqcat("  TxDb: @{object@txdb}\n")
	qqcat("  orgDb: @{object@orgdb}\n")
	qqcat("  Gene sets: @{object@gene_sets_name}\n")

	if(length(object@param)) {
		if(object@param$mode == "basalPlusExt") {
			qqcat("Mode: Basal plus extension\n")
			qqcat("  Proximal: @{object@param$basal_upstream} bp upstream, @{object@param$basal_downstream} bp downstream,\n")
			qqcat("  plus Distal: up to @{object@param$extension} bp\n")

		} else if(object@param$mode == "twoClosest") {
			qqcat("Mode: Two nearest genes\n")
        	qqcat("  within @{object@param$extension} bp\n")
		} else if(object@param$mode == "oneClosest") {
			qqcat("Mode: Single nearest genes\n")
        	qqcat("  within @{object@param$extension} bp\n")
		}
	}
})


# == title
# Get enrichment table
#
# == param
# -object A `GreatObject-class` object returned by `great`.
# -min_hits Minimal number of input regions overlapping to the geneset associated regions.
#
# == details
# Note: adjusted p-values are re-calculated based on ``min_hits``.
#
# == value
# A data frame of enrichment results
#
# == example
# obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
# getEnrichmentTable(obj)
setMethod(f = "getEnrichmentTable",
	signature = "GreatObject",
	definition = function(object, min_hits = 5) {

	tb = object@table
	tb = tb[tb$observed_region_hits >= min_hits, , drop = FALSE]
	tb$p_adjust = p.adjust(tb$p_value, "BH")

	tb
})


# == title
# Get enrichment table
#
# == param
# -object A `GreatObject-class` object returned by `great`.
# -... All passed to `getEnrichmentTable,GreatObject-method`.
#
# == details
# Please use `getEnrichmentTable,GreatObject-method` directly.
#
setMethod(f = "getEnrichmentTables",
	signature = "GreatObject",
	definition = function(object, ...) {
	getEnrichmentTable(object, ...)
})

# == title
# Get region-gene associations
#
# == param
# -object A `GreatObject-class` object returned by `great`.
# -term_id Term ID.
# -use_symbols Whether to use gene symbols
#
# == value
# A `GenomicRanges::GRanges` object. Please the two meta columns are in formats of ``CharacterList``
# and ``IntegerList`` because a region may associate to multiple genes.
#
# == example
# obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
# getRegionGeneAssociations(obj)
setMethod(f = "getRegionGeneAssociations",
	signature = "GreatObject",
	definition = function(object, term_id = NULL, use_symbols = TRUE) {

	gr = object@gr
	extended_tss = object@extended_tss

	gr_mid = gr
	start(gr_mid) = end(gr_mid) = mid(gr)
	ov = findOverlaps(gr_mid, object@background)
	gr = gr[unique(queryHits(ov))]

	if(!is.null(term_id)) {
		if(is.numeric(term_id)) {
			stop_wrap("Do not use numeric index for `term_id`, use the character index.")
		}
		extended_tss = extended_tss[ object@gene_sets[[term_id]] ]
	}

	all_genes = names(extended_tss)

	if(use_symbols) {
		if(grepl("eg.db$", object@orgdb)) {
			map = get_table_from_org_db("egSYMBOL$", object@orgdb)
			map = unlist(as.list(map))
			all_genes2 = map[all_genes]
			l = is.na(all_genes2)
			all_genes2[l] = all_genes[l]
			all_genes = all_genes2
		}
	}

	ov = findOverlaps(gr, extended_tss)

	sr1 = gr[queryHits(ov)]
	sr2 = extended_tss[subjectHits(ov)]

	l_upstream = ( end(sr1) < sr2$tss_position & sr2$tss_strand == "+" ) |
	             ( start(sr1) > sr2$tss_position & sr2$tss_strand == "-" )
	l_downstream = ( end(sr1) < sr2$tss_position & sr2$tss_strand == "-" ) |
	             ( start(sr1) > sr2$tss_position & sr2$tss_strand == "+" )
	dist = pmin( abs(start(sr1) - sr2$tss_position), abs(end(sr1) - sr2$tss_position) )
	dist[l_upstream] = -dist[l_upstream]
	dist[!(l_upstream | l_downstream)] = 0

	# for each region, how many genes it is annotated to
	lt1 = split(subjectHits(ov), queryHits(ov))
	lt_gene = lapply(lt1, function(ind) all_genes[ind])
	gr$annotated_genes = vector("list", length(gr))
	gr$annotated_genes[ as.numeric(names(lt1)) ] = lt_gene

	lt2 = split(dist, queryHits(ov))
	gr$dist_to_TSS = vector("list", length(gr))
	gr$dist_to_TSS[ as.numeric(names(lt2)) ] = lt2

	gr$annotated_genes = CharacterList(gr$annotated_genes)
	gr$dist_to_TSS = IntegerList(gr$dist_to_TSS)

	gr = gr[sapply(gr$annotated_genes, length) > 0]

	return(gr)
})

# == title
# Plot region-gene associations
#
# == param
# -object A `GreatObject-class` object returned by `great`.
# -term_id Term ID.
# -which_plot Which plots to draw? The value should be in ``1, 2, 3``. See "Details" section for explanation.
#
# == details
# There are following figures:  
#
# - Association between regions and genes (``which_plot = 1``).
# - Distribution of distance to TSS (``which_plot = 2``).
# - Distribution of absolute distance to TSS (``which_plot = 3``).
#
# == example
# obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
# plotRegionGeneAssociations(obj)
setMethod(f = "plotRegionGeneAssociations",
    signature = "GreatObject",
    definition = function(object, term_id = NULL, which_plot = 1:3) {

    gr = getRegionGeneAssociations(object, term_id = term_id)

    gr_all = getRegionGeneAssociations(object)

    if(!is.null(term_id)) {
        gr_term = getRegionGeneAssociations(object, term_id = term_id)
    } else {
        gr_term = NULL
    }

    plot_great(gr_all, gr_term, which_plot = which_plot, gr_full_len = length(object@gr), term_id = term_id)

})
