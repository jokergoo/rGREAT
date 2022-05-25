


gr = randomRegions(genome = "hg19")

great(gr, "msigdb:h", "hg19")
great(gr, "msigdb:h", "TxDb.Hsapiens.UCSC.hg19.knownGene")
great(gr, "msigdb:h", "RefSeq:hg19")
great(gr, "msigdb:h", "RefSeqSelect:hg19")
great(gr, "msigdb:h", "GREAT:hg19")
great(gr, "msigdb:h", "Gencode_v19")


for(i in 1:nrow(BIOC_ANNO_PKGS)) {
	genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i]

	if(genome == "araTha") {
		gr = randomRegions(seqlengths = seqlengths(getFromNamespace(BIOC_ANNO_PKGS$txdb[i], ns = BIOC_ANNO_PKGS$txdb[i])))
	} else {
		gr = randomRegions(genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	}

	great(gr, "GO:BP", BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	great(gr, "GO:BP", BIOC_ANNO_PKGS$txdb[i])

	
	if(genome != "araTha") {
		try({
			great(gr, "GO:BP", qq("RefSeqCurated:@{genome}"))
			great(gr, "GO:BP", qq("RefSeqSelect:@{genome}"))
			great(gr, "GO:BP", qq("RefSeq:@{genome}"))
		})
	}
}

for(orgdb in unique(BIOC_ANNO_PKGS$orgdb)) {
	i = which(BIOC_ANNO_PKGS$orgdb == orgdb)[1]

	map = get_table_from_orgdb("CHRLOC$", orgdb)
	if(!is.null(map)) {
		if(genome == "araTha") {
			gr = randomRegions(seqlengths = seqlengths(getFromNamespace(BIOC_ANNO_PKGS$txdb[i], ns = BIOC_ANNO_PKGS$txdb[i])))
		} else {
			gr = randomRegions(genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
		}
		great(gr, "GO:BP", orgdb)
	}
}

for(genome in c("hg19", "hg38", "mm9", "mm10")) {
	gr = randomRegions(genome = genome)

	great(gr, "GO:BP", qq("GREAT:@{genome}"))
}

for(v in paste0("v", 18:19)) {
	gr = randomRegions(genome = "hg19")
	great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("v", 20:40)) {
	gr = randomRegions(genome = "hg38")
	great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("vM", 1)) {
	gr = randomRegions(genome = "mm9")
	great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("vM", 2:25)) {
	gr = randomRegions(genome = "mm10")
	great(gr, "GO:BP", qq("gencode_@{v}"))
}


for(dataset in sample(BIOMART[, 1], 10)) {
	genes = getGenesFromBioMart(dataset)
	sl = tapply(end(genes), seqnames(genes), max)
	sl = structure(as.vector(sl), names = names(sl))
	gr = randomRegions(seqlengths = sl)
	great(gr, "GO:BP", biomart_dataset = dataset)
}


gr = randomRegions(genome = "hg19")
great(gr, "GO:BP", "hg19", background = paste0("chr", 1:22))
great(gr, "GO:BP", "hg19", exclude = c("chrX", "chrY"))
great(gr, "MSigDB:H", "hg19", exclude = getGapFromUCSC("hg19"))


