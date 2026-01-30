
rGREAT_env$NCBIGenomeDownloaded = list()

# == title
# Get genome data from NCBI
#
# == param
# -refseq_assembly_accession The RefSeq accession number for the assembly, such as "GCF_000001405.40" for human.
# -return_granges If the assembly is already on chromosome level, it will directly construct a GRanges object where "chromosomes" are only used and chromosome
#   lengths are corrected fitted in its ``seqlengths``. 
#
# == details
# Only protein coding genes are used.
#
# == value
# If ``return_granges`` is set to ``FALSE``, it returns a list of two data frames:
#
# -genome A data frame of several columns.
# -gene A data frame for genes. The first column contains the RefSeq accession numbers of the corresponding contigs. If the genome is assembled on the chromosome level, the first column corresponds to chromosomes. The contig names can be converted to other names with the information in the ``genome`` data frame.
#
# == example
# if(FALSE) {
# getGenomeDataFromNCBI("GCF_000001405.40", return_granges = TRUE)
# getGenomeDataFromNCBI("GCF_000001405.40")
# }
getGenomeDataFromNCBI = function(refseq_assembly_accession, return_granges = FALSE) {
	if(!grepl("\\.\\d+$", refseq_assembly_accession)) {
		stop("RefSeq accession number for the genome must have version number as the suffix.")
	}
	if(!grepl("^GCF", refseq_assembly_accession)) {
		stop("RefSeq accession number should start with 'GCF' prefix.")
	}

	if(is.null(rGREAT_env$NCBIGenomeDownloaded[[refseq_assembly_accession]])) {

		url = qq("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/@{refseq_assembly_accession}/download?include_annotation_type=GENOME_GFF,SEQUENCE_REPORT&filename=@{refseq_assembly_accession}.zip")

		tmpdir = tempdir()
		destfile = qq("@{tmpdir}/@{refseq_assembly_accession}.zip")
		download.file(url, destfile = destfile, quiet = FALSE)
		on.exit(file.remove(destfile))

		if(Sys.which("unzip") == "") {
			stop("unzip command should be supported.")
		}
		if(Sys.which("perl") == "") {
			stop("Perl should be installed.")
		}

		con = pipe(qq("unzip -p \"@{destfile}\" ncbi_dataset/data/assembly_data_report.jsonl"))
		annotation_info = fromJSON(paste(readLines(con), collapse = "\n"))
		close(con)

		con = pipe(qq("unzip -p \"@{destfile}\" ncbi_dataset/data/@{refseq_assembly_accession}/sequence_report.jsonl"))
		sequence_info = lapply(readLines(con), fromJSON)
		close(con)

		perl_script = system.file("extdata", "download_ncbi_genome_data.pl", package = "rGREAT")
		con = pipe(qq("unzip -p \"@{destfile}\" ncbi_dataset/data/@{refseq_assembly_accession}/genomic.gff | perl \"@{perl_script}\""))
		tb = read.table(con,  sep = "\t")

		colnames(tb) = c("RefSeq Contig Accession", "Start", "End", "Strand", "EntreZ ID")
		tb[, 5] = as.character(tb[, 5])

		rGREAT_env$NCBIGenomeDownloaded[[refseq_assembly_accession]] = list(annotation_info = annotation_info, sequence_info = sequence_info, tb = tb)
	} else {
		annotation_info = rGREAT_env$NCBIGenomeDownloaded[[refseq_assembly_accession]]$annotation_info
		sequence_info = rGREAT_env$NCBIGenomeDownloaded[[refseq_assembly_accession]]$sequence_info
		tb = rGREAT_env$NCBIGenomeDownloaded[[refseq_assembly_accession]]$tb
	}

	cn = c("genbankAccession", "refseqAccession", "chrName", "ucscStyleName", "length")
	if(annotation_info$assemblyInfo$assemblyLevel == "Chromosome") {
		l = sapply(sequence_info, function(x) x$role == "assembled-molecule")
		tb_genome = do.call(rbind, lapply(sequence_info[l], function(x) {
			x2 = rep(list(NA), length(cn))
			names(x2) = cn
			nm = intersect(cn, names(x))
			x2[nm] = x[nm]
			data.frame(x2)
		}))

		tb = tb[tb[, 1] %in% tb_genome$refseqAccession, , drop = FALSE]

	} else {
		tb_genome = do.call(rbind, lapply(sequence_info, function(x) {
			x2 = rep(list(NA), length(cn))
			names(x2) = cn
			nm = intersect(cn, names(x))
			x2[nm] = x[nm]
			data.frame(x2)
		}))
	}

	if(return_granges) {
		if(annotation_info$assemblyInfo$assemblyLevel == "Chromosome") {
			map = structure(tb_genome$chrName, names = tb_genome$refseqAccession)
			tb[, 1] = map[tb[, 1]]

			gene = GRanges(seqnames = tb[, 1], ranges = IRanges(tb[, 2], tb[, 3]), strand = tb[, 4], gene_id = tb[, 5])

			sl = structure(tb_genome$length, names = tb_genome$chrName)
			seqlengths(gene) = sl

			genome(gene) = paste0(annotation_info$organism$organismName, " / ", annotation_info$accession)

			return(gene)
		} else {
			stop_wrap("The genome is not assembled on the chromosome-level. Cannot automatically format it into a GRanges object. Please set `return_granges` to FALSE and construct it manually with the `genome` and `gene` data frames.")
		}
	}

	rownames(tb) = NULL
	lt = list(genome = tb_genome, gene = tb)

	return(lt)
}


