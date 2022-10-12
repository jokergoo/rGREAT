

# == title
# Class for local GREAT analysis
#
# == details
# `great` returns A `GreatObject-class` object.
#
GreatObject = setClass("GreatObject",
    slots = list(table = "data.frame",
    	         gr = "GRanges",
    	         n_total = "numeric",
    	         gene_sets = "list",
    	         extended_tss = "GRanges",
    	         n_gene_gr = "numeric",
    	         gene_sets_name = "character",
    	         background = "GRanges",
    	         background_type = "character",
    	         param = "list")
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


# == title
# Perform GREAT analysis
#
# == param
# -gr A `GenomicRanges::GRanges` object. This is the input regions. It is important to keep consistent for the chromosome names
#     of the input regions and the internal TSS regions. Use `getTSS` to see the format of internal TSS regions.
# -gene_sets A single string of defautly supported gene sets collections (see the full list in "Genesets" section), or a named list of vectors where each vector correspond to a gene set.
# -tss_source Source of TSS. See "TSS" section.
# -biomart_dataset The value should be in ``BioMartGOGeneSets::supportedOrganisms``. 
# -min_gene_set_size Minimal size of gene sets.
# -mode The mode to extend genes. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See `extendTSS` for details.
# -basal_upstream In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.
# -basal_downstream In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.
# -extension Extensions from the basal domains.
# -extended_tss If your organism is not defaultly supported, you can first prepare one by `extendTSSFromDataFrame` or `extendTSS`,
#      and set the object to this argument. Please see more examples in the vignette.
# -background Background regions. The value can also be a vector of chromosome names.
# -exclude Regions that are excluded from analysis such as gap regions (which can be get by `getGapFromUCSC`). The value can also be a vector of chromosome names.
#     It also allows a special character value ``"gap"`` so that gap regions for corresponding organism will be removed from the analysis.
# -cores Number of cores to use.
# -verbose Whether to print messages.
#
# == details
# When ``background`` or ``exclude`` is set, the analysis is restricted in the background regions, still by using Binomial method. Note
# this is different from the original GREAT method which uses Fisher's exact test if background regions is set. See `submitGreatJob` for explanations.
#
# By default, gap regions are excluded from the analysis.
#
# == TSS
# rGREAT supports TSS from many organisms. The value of ``tss_source`` should be encoded in a special format:
#
# - Name of ``TxDb.*`` packages. Supported packages are in ``rGREAT:::BIOC_ANNO_PKGS$txdb``.
# - Genome version of the organism, e.g. "hg19". Then the corresponding TxDb will be used.
# - In a format of ``RefSeqCurated:$genome`` where ``$genome`` is the genome version of an organism, such as hg19. RefSeqCurated subset will be used.
# - In a format of ``RefSeqSelect:$genome`` where ``$genome`` is the genome version of an organism, such as hg19. RefSeqSelect subset will be used.
# - In a format of ``Gencode_v$version`` where ``$version`` is gencode version, such as 19 (for human) or M21 for mouse. Gencode protein coding genes will be used.
# - In a format of ``GREAT:$genome``, where ``$genome`` can only be mm9, mm10, hg19, hg38. The TSS from GREAT will be used.
#
# == Genesets
#
# rGREAT supports the following built-in GO gene sets for all organisms (note "GO:" can be omitted):
#
# -"GO:BP": Biological Process, from GO.db package.
# -"GO:CC": Cellular Component, from GO.db package.
# -"GO:MP": Molecular Function, from GO.db pacakge.
#
# rGREAT also supports built-in gene sets collections from MSigDB (note this is only for human, "msigdb:" can be omitted):
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
# If the defaultly supported TxDb is used, Entrez gene ID is always used as the main gene ID. If you provide a self-defined
# ``gene_sets`` or ``extended_tss``, you need to make sure they two have the same gene ID types.
#
# == BioMart
#
# rGREAT supports a large number of organisms of which the information is retrieved from Ensembl BioMart. The name of a BioMart dataset
# can be assigned to argument ``biomart_dataset``. All supported organisms can be found with ``BioMartGOGeneSets::supportedOrganisms``.
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
# gr = randomRegions(genome = "hg19")
# res = great(gr, "MSigDB:H", "txdb:hg19")
# res = great(gr, "MSigDB:H", "TxDb.Hsapiens.UCSC.hg19.knownGene")
# res = great(gr, "MSigDB:H", "RefSeq:hg19")
# res = great(gr, "MSigDB:H", "GREAT:hg19")
# res = great(gr, "MSigDB:H", "Gencode_v19")
# res = great(gr, "GO:BP", "hsapiens_gene_ensembl")
# }
great = function(gr, gene_sets, tss_source, biomart_dataset = NULL,
	min_gene_set_size = 5, mode = "basalPlusExt", basal_upstream = 5000, 
	basal_downstream = 1000, extension = 1000000,
	extended_tss = NULL, background = NULL, exclude = "gap",
	cores = 1, verbose = great_opt$verbose) {

	param = list()
	gene_sets_name = ""
	if(is.null(extended_tss)) {
		if(!is.null(biomart_dataset)) {
			
			if(is.null(rGREAT_env$BIOMART)) {
				check_pkg("BioMartGOGeneSets", bioc = TRUE, github = "jokergoo")
				rGREAT_env$BIOMART = BioMartGOGeneSets::supportedOrganisms(html = FALSE)
			}
			if(is.null(rGREAT_env$BIOMART$genome)) {
				rGREAT_env$BIOMART$genome = BIOMART_GENOME[rGREAT_env$BIOMART[, 1]]
				rGREAT_env$BIOMART$genome[is.na(rGREAT_env$BIOMART$genome)] = ""
			}

			BIOMART = rGREAT_env$BIOMART

			biomart_dataset = tolower(biomart_dataset)
			if(!biomart_dataset %in% BIOMART[, 1]) {
				i = grep(paste0("\\b", biomart_dataset, "\\b"), BIOMART[, 1], ignore.case = TRUE)
				if(length(i) == 1) {
					biomart_dataset = BIOMART[i, 1]
				} else {
					i = grep(paste0("\\b", biomart_dataset, "\\b"), BIOMART[, "genome"], ignore.case = TRUE)
					if(length(i) == 1) {
						biomart_dataset = BIOMART[i, 1]
					} else {
						i = grep(paste0("\\b", biomart_dataset, "\\b"), BIOMART[, "description"], ignore.case = TRUE)
						if(length(i) == 1) {
							biomart_dataset = BIOMART[i, 1]
						} else {
							stop_wrap(qq("Cannot find biomart dataset: @{biomart_dataset}."))
						}
					}
				}
			}

			if(verbose) {
				message(qq("* get gene and GO genesets from biomart (dataset: @{biomart_dataset})."))
			}

			gene_sets = tolower(gene_sets)
			if(gene_sets %in% c("bp", "go:bp")) {
				gene_sets_name = "GO:BP"
				gene_sets = getGeneSetsFromBioMart(biomart_dataset, "bp")
			} else if(gene_sets %in% c("cc", "go:cc")) {
				gene_sets_name = "GO:CC"
				gene_sets = getGeneSetsFromBioMart(biomart_dataset, "cc")
			} else if(gene_sets %in% c("mf", "go:mf")) {
				gene_sets_name = "GO:MF"
				gene_sets = getGeneSetsFromBioMart(biomart_dataset, "mf")
			} else {
				stop_wrap("When `biomart_dataset` is set, `gene_sets` can only be one of 'GO:BP/GO:CC/GO:MF'.")
			}

			if(great_opt$test) {
				gene_sets = gene_sets[order(-sapply(gene_sets, length))[1:10]]
			}

			i = which(BIOMART[, 1] == biomart_dataset)
			genome = BIOMART[i, "genome"]

			genes = getGenesFromBioMart(biomart_dataset, filter = TRUE)

			if(genome == "") {
				sl = tapply(end(genes), seqnames(genes), max)
				sl = structure(as.vector(sl), names = names(sl))
			} else {
				allg = unique(unlist(gene_sets))
				all_seqlevels = unique(seqnames(genes)[genes$gene_id %in% allg])

				oe = try(sl <- getChromInfoFromUCSC(genome, seqlevels = all_seqlevels), silent = TRUE)
				
				if(inherits(oe, "try-error")) {
					sl = tapply(end(genes), seqnames(genes), max)
					sl = structure(as.vector(sl), names = names(sl))
				} else {
					names(sl) = gsub("^chr", "", names(sl))
				}
				# there might be inconsistency between UCSC and Ensembl...
				if(length(intersect(seqlevels(genes), names(sl))) == 0) {
					sl = tapply(end(genes), seqnames(genes), max)
					sl = structure(as.vector(sl), names = names(sl))
				}
			}

			if(verbose) {
				message_str = qq("* background chromosomes from BioMart are '@{paste(names(sl), collapse=', ')}'.")
				message_str = strwrap(message_str, width = 90)
				if(length(message_str) > 1) {
					message_str[2:length(message_str)]= paste0("    ", message_str[2:length(message_str)])
				}
				message(paste(message_str, collapse = "\n"))
			}

			extended_tss = extendTSS(genes, seqlengths = sl, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)

			param$genome = BIOMART[i, "description"]
			if(genome != "") {
				param$genome = paste0(param$genome, '; ', genome)
			}
			
			param$tss_source = biomart_dataset

			gr = GRanges(seqnames = gsub("^chr", "", as.vector(seqnames(gr))), ranges = ranges(gr))
			if(inherits(background, "GRanges")) {
				background = GRanges(seqnames = gsub("^chr", "", as.vector(seqnames(background))), ranges = ranges(background))
			}
			if(inherits(exclude, "GRanges")) {
				exclude = GRanges(seqnames = gsub("^chr", "", as.vector(seqnames(exclude))), ranges = ranges(exclude))
			} else if(identical(exclude, "gap") || identical(exclude, "gaps")) {
				oe = try(exclude <- getGapFromUCSC(genome), silent = TRUE)
				if(inherits(oe, "try-error")) {
					exclude = NULL
				} else {
					new_seqnames = gsub("^chr", "", as.vector(seqnames(exclude)))
					l = new_seqnames %in% seqlevels(extended_tss)
					exclude = GRanges(seqnames = new_seqnames[l], ranges = ranges(exclude)[l])
					attr(exclude, "exclude_is_gap") = TRUE
				}
			}
		} else {

			tss_source = parse_tss_source(tss_source)

			orgdb = get_orgdb_from_genome_version(tss_source$genome)
			if(is.character(gene_sets) && is.atomic(gene_sets)) {
				if(orgdb == "") {
					stop_wrap("There is no pre-defined gene sets for this organism.")
				}
			}

			if(tss_source$category == "TxDb") {
				if(verbose) {
					message("* TSS source: TxDb.")
				}
				extended_tss = extendTSSFromTxDb(tss_source$source, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = tss_source$source 
			} else if(tss_source$category == "OrgDb") {
				if(verbose) {
					message("* TSS source: OrgDb")
				}
				extended_tss = extendTSSFromOrgDb(tss_source$source, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = tss_source$source 
			} else if(tss_source$category == "Gencode") {
				if(verbose) {
					message("* TSS source: Gencode.")
				}
				tss = getGenesFromGencode(tss_source$source)
				extended_tss = extendTSS(tss, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = paste0("Gencode_", tss_source$source)
			} else if(tss_source$category == "RefSeq") {
				if(verbose) {
					message("* TSS source: RefSeqSelect.")
				}
				tss = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqSelect")
				extended_tss = extendTSS(tss, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = "RefSeq"
			} else if(tss_source$category == "RefSeqCurated") {
				if(verbose) {
					message("* TSS source: RefSeqCurated.")
				}
				tss = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqCurated")
				extended_tss = extendTSS(tss, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = "RefSeqCurated"
			} else if(tss_source$category == "RefSeqSelect") {
				if(verbose) {
					message("* TSS source: RefSeqSelect.")
				}
				tss = getRefSeqGenesFromUCSC(tss_source$genome, subset = "RefSeqSelect")
				extended_tss = extendTSS(tss, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = "RefSeqSelect"
			} else if(tss_source$category == "GREAT") {
				if(verbose) {
					message("* TSS source: GREAT.")
				}
				tss = getGREATDefaultTSS(tss_source$genome)
				extended_tss = extendTSS(tss, mode = mode, basal_upstream = basal_upstream,
					basal_downstream = basal_downstream, extension = extension)
				param$tss_source = "GREAT"
			} else {
				stop_wrap("Wrong `tss_source`.")
			}

			param$genome = tss_source$genome
		}

		seqlevelsStyle(extended_tss) = seqlevelsStyle(gr)[1]

		# check genomic end coordinate and chr length
		sl = seqlengths(extended_tss)
		if(!is.na(sl[1])) {
			gre = tapply(end(gr), seqnames(gr), max)
			gre[is.na(gre)] = 0

			for(cn in intersect(names(sl), names(gre))) {
				if(sl[cn] < gre[cn]) {
					warning_wrap(qq("Inconsistent coordinate on '@{cn}': Max end position (@{gre[cn]}) is larger than chromosome length (@{sl[cn]})."))
				}
			}
		}
		
	} else {
		mt = attributes(metadata(extended_tss))
		
		if(is.null(mt$genome)) {
			param$genome = "unknown"
		} else {
			param$genome = mt$genome
		}
		param$tss_source = "self-provided"

		seqlevelsStyle(extended_tss) = seqlevelsStyle(gr)[1]
	}

	param$orgdb = get_orgdb_from_genome_version(param$genome)

	extended_tss = sort(extended_tss)
	mode_param = attr(extended_tss, "mode")
	param = c(param, mode_param)

	mt = attributes(metadata(extended_tss))
	if(is.character(gene_sets) && is.atomic(gene_sets)) {

		gene_sets_name = gene_sets
		mt = attributes(metadata(extended_tss))

		gene_sets = get_defaultly_suppported_gene_sets(gene_sets, mt, verbose = verbose,
			gene_id = extended_tss$gene_id)

		if(is.null(gene_sets)) {
			stop_wrap("Cannot find the pre-defined gene sets for the specified organism. Please set `gene_sets` as a named list of vectors where each vector corresponds to a gene set. The gene ID type should be the same as in `txdb` or `extended_tss`. If the defaultly supported TxDb package is used, Entrez gene ID is the main gene ID.")
		}
		
	} else {
		if(gene_sets_name == "") {
			gene_sets_name = "self-provided"
		}
	}

	gene_sets = lapply(gene_sets, unique)

	if(!"gene_id" %in% colnames(mcols(extended_tss))) {
		stop_wrap("`extended_tss` must have a meta column `gene_id`")
	}
	names(extended_tss) = extended_tss$gene_id

	if(verbose) {
		message("* check gene ID type in `gene_sets` and in `extended_tss`.")
	}

	v1 = unique(unlist(gene_sets))
	if(length(intersect(v1, extended_tss$gene_id))/length(v1) <= 0) {
		stop_wrap(qq("It seems the gene ID type in `gene_sets` (e.g. '@{gene_sets[[1]][1]}') is different from in `extended_tss` (e.g. '@{names(extended_tss)[1]}')."))
	}
	
	sl = seqlengths(extended_tss)
	if(any(is.na(sl))) {
		stop_wrap("`extended_tss` must have `seqlengths` been set.")
	}

	gr_genome = GRanges(seqnames = names(sl), ranges = IRanges(1, sl))
	background_type = "whole genome"
	if(is.null(background)) {
		background = gr_genome

		if(verbose) {
			message("* use whole genome as background.")
		}
	} else {
		if(verbose) {
			message("* use self-defined background regions.")
		}
		if(inherits(background, "character")) {
			background_chr = background
			if(!is.null(biomart_dataset)) {
				background_chr = gsub("^chr", "", background_chr)
			}
			background = gr_genome[seqnames(gr_genome) %in% background_chr]
		}
		background_type = "self-provided"
	}
	strand(background) = "*"

	if(!is.null(exclude)) {
		if(identical(exclude, "gap") || identical(exclude, "gaps")) {
			
			if(!is.null(param$genome)) {
				genome = param$genome
				oe = try(exclude <- getGapFromUCSC(genome, seqnames = seqlevels(background)), silent = TRUE)
				if(!inherits(oe, "try-error")) {
					strand(exclude) = "*"
					background = setdiff(background, exclude)

					if(background_type == "whole genome") {
						background_type = "whole genome excluding gaps"
					}
				}
			}

		} else if(inherits(exclude, "character")) {
			background = background[!seqnames(background) %in% exclude]
			background_type = "self-provided"
		} else {
			strand(exclude) = "*"
			background = setdiff(background, exclude)

			if(background_type == "whole genome") {
				if(!is.null(attr(exclude, "exclude_is_gap"))) {
					background_type = "whole genome excluding gaps"
				} else {
					background_type = "self-provided"
				}
			} else {
				background_type = "self-provided"
			}
		}

		if(verbose){
			message("* remove excluded regions from background.")
		}
	}

	background = reduce(background)
	background_total_length = sum(width(background))

	strand(gr) = "*"

	gr_origin = gr

	gr_mid = gr
	start(gr_mid) = end(gr_mid) = mid(gr)
	gr = gr_mid

	if(verbose) {
		message("* overlap `gr` to background regions (based on midpoint).")
	}

	gr = intersect(gr, background)  # possible seqlevels are different due to different versions of organisms

	n_total = length(gr)

	if(verbose) {
		message(qq("* in total @{n_total} `gr`."))
	}

	if(verbose) {
		message(qq("* overlap extended TSS to background regions."))
	}

	ov2 = findOverlaps(extended_tss, background)
	r1 = extended_tss[queryHits(ov2)]
	r2 = background[subjectHits(ov2)]

	# note in extended_tss, there might be same gene in multi rows, due to a gene is segmented by background
	extended_tss2 = pintersect(r1, r2)
	extended_tss2$gene_id = r1$gene_id
	names(extended_tss2) = r1$gene_id
	extended_tss = extended_tss2
	extended_tss$hit = NULL

	# note, after overlap to background, an extended TSS may be split into several regions
	if(verbose) {
		message("* check which genes are in the gene sets.")
	}
	gene_sets = lapply(gene_sets, function(x) intersect(x, names(extended_tss)))
	
	if(verbose) {
		message(qq("* only take gene sets with size >= @{min_gene_set_size}."))
	}
	l = sapply(gene_sets, length) >= min_gene_set_size
	gene_sets = gene_sets[l]
	gene_sets_ind = lapply(gene_sets, function(x) which(names(extended_tss) %in% x))
	
	n_set = length(gene_sets)

	if(n_set == 0) {
		stop_wrap("No gene set left.")
	}
	
	if(verbose) {
		message(qq("* in total @{n_set} gene sets."))
	}

	p = numeric(n_set)
	n_obs = numeric(n_set)
	n_exp = numeric(n_set)
	fold = numeric(n_set)
	prop = numeric(n_set)
	mean_tss_dist = numeric(n_set)

	gene_hits = numeric(n_set)

	if(verbose) {
		message(qq("* overlap `gr` to every extended TSS."))
	}
	ov = findOverlaps(gr, extended_tss)

	if(verbose) {
		message("* perform binomial test for each biological term.")
	}
	pb = progress::progress_bar$new(total = n_set)

	all_chr = as.vector(seqnames(extended_tss))
	all_s = start(extended_tss)
	all_e = end(extended_tss)
	tss_position = extended_tss$tss_position
	gr_start = start(gr)

	registerDoParallel(cores)
	lt <- foreach(i = seq_len(n_set)) %dopar% {
		if(verbose && cores <= 1) {
			pb$tick()
		}
		ind = gene_sets_ind[[i]]
		
		if(length(ind) == 0) {
			width_fgr = 0
		} else {
			sl = split(all_s[ind], all_chr[ind])
			el = split(all_e[ind], all_chr[ind])
			width_fgr = sum(sapply(seq_along(sl), function(i) {
				reduce_by_start_and_end(sl[[i]], el[[i]])
			}))
		}
		prop = width_fgr/background_total_length
	
		l = subjectHits(ov) %in% ind
		n_hits = length(unique(queryHits(ov)[l]))

		mean_tss_dist = mean(abs(tss_position[subjectHits(ov)[l]] - gr_start[queryHits(ov)[l]]))

		if(n_hits == 0) {
			p = 1
		} else {
			p = 1 - pbinom(n_hits - 1, n_total, prop)
		}
		n_obs = n_hits
		n_exp = prop*n_total

		gene_hits = length(unique(extended_tss$gene_id[subjectHits(ov)[l]]))

		return(list(prop = prop, n_obs = n_obs, n_exp = n_exp, p = p, mean_tss_dist = mean_tss_dist, gene_hits = gene_hits))
	}
	stopImplicitCluster()

	prop = sapply(lt, function(x) x$prop)
	n_obs = sapply(lt, function(x) x$n_obs)
	n_exp = sapply(lt, function(x) x$n_exp)
	p = sapply(lt, function(x) x$p)
	mean_tss_dist = sapply(lt, function(x) x$mean_tss_dist)
	gene_hits = sapply(lt, function(x) x$gene_hits)

	fold = n_obs/n_exp

	df = data.frame(id = names(gene_sets),
			   genome_fraction = prop,
		       observed_region_hits = n_obs,
		       fold_enrichment = fold,
		       p_value = p,
		       p_adjust = p.adjust(p, "BH"),
		       mean_tss_dist = round(mean_tss_dist),
		       observed_gene_hits = gene_hits,
		       gene_set_size = sapply(gene_sets, length))

	n_gene_total = length(unique(extended_tss$gene_id))
	n_gene_gr = length(unique(extended_tss[subjectHits(ov)]$gene_id))
	df$fold_enrichment_hyper = df$observed_gene_hits/(df$gene_set_size*n_gene_gr/n_gene_total)
	df$p_value_hyper = 1 - phyper(df$observed_gene_hits - 1, df$gene_set_size, n_gene_total - df$gene_set_size, n_gene_gr)
	df$p_adjust_hyper = p.adjust(df$p_value_hyper, "BH")

	df = df[order(df$p_adjust, df$p_value, -df$fold_enrichment), , drop = FALSE]

	if(all(grepl("^GO:\\d+$", df$id))) {
		check_pkg("GO.db", bioc = TRUE)
		terms = Term(GO.db::GOTERM)

		desc = unname(terms[df$id])
		desc[is.na(desc)] = ""
		df = cbind(data.frame(id = df$id, description = desc),
			       df[, -1])
	}
	rownames(df) = NULL

	obj = GreatObject(
		table = df,
		gr = gr_origin,
		n_total = n_total,
		gene_sets = gene_sets,
		gene_sets_name = gene_sets_name,
		extended_tss = extended_tss,  # has been intersected with background,
		n_gene_gr = n_gene_gr,
		background = background,
		background_type = background_type, 
		param = param
	)

	return(obj)
}


# hg19, txdb:hg19, TxDb.Hsapiens.UCSC.hg19.knownGene, great:hg19, refseq:hg19, refseqcurated:hg19, gencode_v19
# human, txdb:chicken, refseq:mouse
parse_tss_source = function(x) {
	if(grepl("^TxDb\\.", x)) {
		i = which(BIOC_ANNO_PKGS$txdb == x)
		if(length(i) == 0) {
			stop_wrap("TxDb: '@{x}' is not supported.")
		}
		list(category = "TxDb", source = x, genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	} else if(grepl("^txdb:", x, ignore.case = TRUE)) {
		i = detect_txdb(gsub("^txdb:", "", x, ignore.case = TRUE))
		if(length(i) == 0) {
			stop_wrap("Cannot parse `tss_source`.")
		}
		i = i[1]
		list(category = "TxDb", source = BIOC_ANNO_PKGS$txdb[i], genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	} else if(grepl("^org\\.", x)) {
		i = which(BIOC_ANNO_PKGS$orgdb == x)
		if(length(i) == 0) {
			stop_wrap("OrgDb: '@{x}' is not supported.")
		}
		i = i[1]
		list(category = "OrgDb", source = BIOC_ANNO_PKGS$orgdb[i], genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	} else if(grepl("^orgdb:", x, ignore.case = TRUE)) {
		i = detect_txdb(gsub("^orgdb:", "", x, ignore.case = TRUE))
		if(length(i) == 0) {
			stop_wrap("Cannot parse `tss_source`.")
		}
		i = i[1]
		list(category = "OrgDb", source = BIOC_ANNO_PKGS$orgdb[i], genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	} else if(grepl("^great:", x, ignore.case = TRUE)) {
		list(category = "GREAT", genome = gsub("^great:", "", x, ignore.case = TRUE))
	} else if(grepl("^refseq:", x, ignore.case = TRUE)) {
		list(category = "RefSeq", genome = gsub("^refseq:", "", x, ignore.case = TRUE))
	} else if(grepl("^refseqcurated:", x, ignore.case = TRUE)) {
		list(category = "RefSeqCurated", genome = gsub("^refseqcurated:", "", x, ignore.case = TRUE))
	} else if(grepl("^refseqselect:", x, ignore.case = TRUE)) {
		list(category = "RefSeqSelect", genome = gsub("^refseqselect:", "", x, ignore.case = TRUE))
	} else if(grepl("^gencode", x, ignore.case = TRUE)) {
		version = gsub("^gencode_?", "", x, ignore.case = TRUE)
		version = gsub("^V", "v", version)
		if(!grepl("^v", version)) {
			version = paste0("v", version)
		}
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
		list(category = "Gencode", source = version, genome = genome)
	} else if(grepl("^[a-zA-Z]+\\d*$", x)) {
		i = detect_txdb(x)
		if(length(i) == 0) {
			stop_wrap("Cannot parse `tss_source`.")
		}
		i = i[1]
		list(category = "TxDb", source = BIOC_ANNO_PKGS$txdb[i], genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	} else {
		stop_wrap("Cannot parse `tss_source`.")
	}
}


get_orgdb_from_genome_version = function(x) {
	if(x == "unknown") {
		""
	} else {
		map = tapply(BIOC_ANNO_PKGS$orgdb, BIOC_ANNO_PKGS$genome_version_in_txdb, unique)
		v = map[x]
		if(is.na(v)) {
			""
		} else {
			v
		}
	}
}

# support two types
# -GO (GO:BP, BP, GO:CC, CC, GO:MF, MF)
# -MsigDB
get_defaultly_suppported_gene_sets = function(name, mt, verbose = great_opt$verbose, gene_id = NULL) {
	if(length(mt) == 0) {
		return(NULL)
	}
	if(is.null(mt$genome)) {
		return(NULL)
	}

	gene_id_type = mt$gene_id_type
	if(is.null(gene_id_type)) {
		if(!is.null(gene_id)) {
			gene_id_type = guess_gene_id_type(gene_id)
		}
	}

	orgdb = BIOC_ANNO_PKGS$orgdb[ BIOC_ANNO_PKGS$genome_version_in_txdb == mt$genome ][1]
	if(is.na(orgdb)) {
		return(NULL)
	}

	name = tolower(name)

	if(grepl("^c\\d", name) || name == "h") {
		name = paste0("msigdb:", name)
	}

	if(name %in% c("go:bp", "bp", "go:cc", "cc", "go:mf", "mf")) {
		lt = as.list(get_table_from_orgdb("GO2ALL.*S$", orgdb))

		check_pkg("GO.db", bioc = TRUE)

		ontology = Ontology(GO.db::GOTERM)

		if(name %in% c("go:bp", "bp")) {
			l = ontology[names(lt)] == "BP"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:BP ontology with @{length(gene_sets)} gene sets (source: @{orgdb})."))
			}
		} else if(name %in% c("go:mf", "mf")) {
			l = ontology[names(lt)] == "MF"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:MF ontology with @{length(gene_sets)} gene sets (source: @{orgdb})."))
			}
		} else if(name %in% c("go:cc", "cc")) {
			l = ontology[names(lt)] == "CC"
			l[is.na(l)] = FALSE
			gene_sets = lt[l]

			if(verbose) {
				message(qq("* use GO:CC ontology with @{length(gene_sets)} gene sets (source: @{orgdb})."))
			}
		}

		if(great_opt$test) {
			gene_sets = gene_sets[order(-sapply(gene_sets, length))[1:10]]
		}

		if(!is.null(mt$orgdb)) {
			return(gene_sets)
		}

		if(is.null(gene_id_type)) {
			return(gene_sets)
		} else if(gene_id_type == "ENTREZ" || gene_id_type == "Entrez Gene ID" || gene_id_type == "SGD Gene ID" || gene_id_type == "TAIR ID") {
			return(gene_sets)
		} else {
			if(gene_id_type == "ENSEMBL" || gene_id_type == "Ensembl gene ID") {
				map = get_table_from_orgdb("ENSEMBL$", orgdb)
			} else if(gene_id_type == "SYMBOL") {
				map = get_table_from_orgdb("SYMBOL$", orgdb)
			} else if(gene_id_type == "REFSEQ") {
				map = get_table_from_orgdb("REFSEQ$", orgdb)
			} else {
				stop_wrap("Wrong gene_id_type.")
			}

			map = unlist(as.list(map))
			gene_sets = lapply(gene_sets, function(x) {
				x = map[x]
				unique(x[!is.na(x)])
			})
			gene_sets = gene_sets[sapply(gene_sets, length) > 0]
			return(gene_sets)
		}

	} else if(grepl("msigdb", name)) {
		if(!grepl("^hg\\d", mt$genome)) {
			stop_wrap("MSigDB only supports human.")
		}

		name = gsub(":", ".", name)
		name = gsub("^msigdb\\.?", "", name)
		if(name %in% names(MSIGDB)) {
			gene_sets = readRDS(get_url(MSIGDB[[name]]))
		} else {
			stop_wrap("Wrong MSigDB gene set collection. All supported gene sets are in `names(rGREAT:::MSIGDB)`.")
		}

		orgdb = "org.Hs.eg.db"
		if(is.null(gene_id_type)) {
			return(gene_sets)
		} else if(gene_id_type == "ENTREZ" || gene_id_type == "Entrez Gene ID" || gene_id_type == "SGD Gene ID" || gene_id_type == "TAIR ID") {
			return(gene_sets)
		} else {
			if(gene_id_type == "ENSEMBL" || gene_id_type == "Ensembl gene ID") {
				map = get_table_from_orgdb("ENSEMBL$", orgdb)
			} else if(gene_id_type == "SYMBOL") {
				map = get_table_from_orgdb("SYMBOL$", orgdb)
			} else if(gene_id_type == "REFSEQ") {
				map = get_table_from_orgdb("REFSEQ$", orgdb)
			} else {
				stop_wrap("Wrong gene_id_type.")
			}

			map = unlist(as.list(map))
			gene_sets = lapply(gene_sets, function(x) {
				x = map[x]
				unique(x[!is.na(x)])
			})
			gene_sets = gene_sets[sapply(gene_sets, length) > 0]
			return(gene_sets)
		}

		if(verbose) {
			message(qq("* use @{name} collection with @{length(gene_sets)} gene sets (source: MSigDb)."))
		}
		return(gene_sets)
	} else {
		stop_wrap(qq("Gene sets '@{name}' is not supported."))
	}
}

MSIGDB = c(
	"c1" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c1.all.v7.5.1.entrez.gmt.rds",
	"c2" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.all.v7.5.1.entrez.gmt.rds",
	"c2.cgp" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cgp.v7.5.1.entrez.gmt.rds",
	"c2.cp.kegg" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cp.kegg.v7.5.1.entrez.gmt.rds",
	"c2.cp.pid" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cp.pid.v7.5.1.entrez.gmt.rds",
	"c2.cp.reactome" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cp.reactome.v7.5.1.entrez.gmt.rds",
	"c2.cp" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cp.v7.5.1.entrez.gmt.rds",
	"c2.cp.wikipathways" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c2.cp.wikipathways.v7.5.1.entrez.gmt.rds",
	"c3" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.all.v7.5.1.entrez.gmt.rds",
	"c3.mir.mir_legacy" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.mir.mir_legacy.v7.5.1.entrez.gmt.rds",
	"c3.mir.mirdb" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.mir.mirdb.v7.5.1.entrez.gmt.rds",
	"c3.mir" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.mir.v7.5.1.entrez.gmt.rds",
	"c3.tft.gtrd" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.tft.gtrd.v7.5.1.entrez.gmt.rds",
	"c3.tft.tft_legacy" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.tft.tft_legacy.v7.5.1.entrez.gmt.rds",
	"c3.tft" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c3.tft.v7.5.1.entrez.gmt.rds",
	"c4" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c4.all.v7.5.1.entrez.gmt.rds",
	"c4.cgn" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c4.cgn.v7.5.1.entrez.gmt.rds",
	"c4.cm" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c4.cm.v7.5.1.entrez.gmt.rds",
	"c5" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.all.v7.5.1.entrez.gmt.rds",
	"c5.go.bp" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.go.bp.v7.5.1.entrez.gmt.rds",
	"c5.go.cc" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.go.cc.v7.5.1.entrez.gmt.rds",
	"c5.go.mf" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.go.mf.v7.5.1.entrez.gmt.rds",
	"c5.go" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.go.v7.5.1.entrez.gmt.rds",
	"c5.hpo" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c5.hpo.v7.5.1.entrez.gmt.rds",
	"c6" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c6.all.v7.5.1.entrez.gmt.rds",
	"c7" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c7.all.v7.5.1.entrez.gmt.rds",
	"c7.immunesigdb" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c7.immunesigdb.v7.5.1.entrez.gmt.rds",
	"c7.vax" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c7.vax.v7.5.1.entrez.gmt.rds",
	"c8" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/c8.all.v7.5.1.entrez.gmt.rds",
	"h" = "https://jokergoo.github.io/rGREAT_genesets/msigdb/h.all.v7.5.1.entrez.gmt.rds"
)


setMethod(f = "show",
	signature = "GreatObject",
	definition = function(object) {

	qqcat("@{object@n_total} regions are associated to @{object@n_gene_gr} genes' extended TSSs.\n")
	qqcat("  TSS source: @{object@param$tss_source}\n")
	qqcat("  Genome: @{object@param$genome}\n")
	if(object@param$orgdb != "") {
		qqcat("  OrgDb: @{object@param$orgdb}\n")
	}
	qqcat("  Gene sets: @{object@gene_sets_name}\n")
	qqcat("  Background: @{object@background_type}\n")

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
# -min_region_hits Minimal number of input regions overlapping to the geneset associated regions.
#
# == details
# Note: adjusted p-values are re-calculated based on ``min_region_hits``.
#
# == value
# A data frame of enrichment results
#
# == example
# obj = readRDS(system.file("extdata", "GreatObject.rds", package = "rGREAT"))
# getEnrichmentTable(obj)
setMethod(f = "getEnrichmentTable",
	signature = "GreatObject",
	definition = function(object, min_region_hits = 5) {

	tb = object@table
	tb = tb[tb$observed_region_hits >= min_region_hits, , drop = FALSE]
	tb$p_adjust = p.adjust(tb$p_value, "BH")
	tb$p_adjust_hyper = p.adjust(tb$p_value_hyper, "BH")

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
# == value
# A data frame of enrichment results
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
# -by_middle_points Whether the distances are calculated from the middle points of input regions?
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
	definition = function(object, term_id = NULL, by_middle_points = FALSE, 
		use_symbols = TRUE) {

	gr = object@gr
	extended_tss = object@extended_tss

	gr_mid = gr
	start(gr_mid) = end(gr_mid) = mid(gr)
	if(by_middle_points) {
		gr = gr_mid
	}
	ov = findOverlaps(gr_mid, object@background)
	gr = gr[unique(queryHits(ov))]

	if(!is.null(term_id)) {
		if(is.numeric(term_id)) {
			stop_wrap("Do not use numeric index for `term_id`, use the character index.")
		}
		extended_tss = extended_tss[ extended_tss$gene_id %in% object@gene_sets[[term_id]] ]
	}

	all_genes = names(extended_tss)

	if(use_symbols) {
		if(grepl("eg.db$", object@param$orgdb)) {
			map = get_table_from_orgdb("egSYMBOL$", object@param$orgdb)
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

    plot_great(gr_all, gr_term, which_plot = which_plot, gr_full_len = object@n_total, term_id = term_id)

})
