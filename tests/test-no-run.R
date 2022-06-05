

test_great = function(...) {
	
	cat("======= test start =========\n")
	res = great(...)
	print(res)
	cat("------- top terms --------\n")
	tb = getEnrichmentTable(res)
	print(head(tb))
}

#### msigdb #####
gr = randomRegions(genome = "hg19")

test_great(gr, "msigdb:h", "hg19")
test_great(gr, "msigdb:h", "TxDb.Hsapiens.UCSC.hg19.knownGene")
test_great(gr, "msigdb:h", "RefSeq:hg19")
test_great(gr, "msigdb:h", "RefSeqSelect:hg19")
test_great(gr, "msigdb:h", "GREAT:hg19")
test_great(gr, "msigdb:h", "Gencode_v19")



great_opt = setGlobalOptions(
	verbose = TRUE,
	test = TRUE
)

###### txdb and refseq ###

for(i in 1:nrow(BIOC_ANNO_PKGS)) {
	genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i]

	if(genome == "araTha") {
		gr = randomRegions(seqlengths = seqlengths(getFromNamespace(BIOC_ANNO_PKGS$txdb[i], ns = BIOC_ANNO_PKGS$txdb[i])))
	} else {
		gr = randomRegions(genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	}

	test_great(gr, "GO:BP", BIOC_ANNO_PKGS$genome_version_in_txdb[i])
	test_great(gr, "GO:BP", BIOC_ANNO_PKGS$txdb[i])

	
	if(genome != "araTha") {
		try({
			test_great(gr, "GO:BP", qq("RefSeqCurated:@{genome}"))
			test_great(gr, "GO:BP", qq("RefSeqSelect:@{genome}"))
			test_great(gr, "GO:BP", qq("RefSeq:@{genome}"))
		})
	}
}


#### orgdb #####

for(orgdb in unique(BIOC_ANNO_PKGS$orgdb)) {
	i = which(BIOC_ANNO_PKGS$orgdb == orgdb)[1]

	map = get_table_from_orgdb("CHRLOC$", orgdb)
	if(!is.null(map)) {
		if(genome == "araTha") {
			gr = randomRegions(seqlengths = seqlengths(getFromNamespace(BIOC_ANNO_PKGS$txdb[i], ns = BIOC_ANNO_PKGS$txdb[i])))
		} else {
			gr = randomRegions(genome = BIOC_ANNO_PKGS$genome_version_in_txdb[i])
		}
		test_great(gr, "GO:BP", orgdb)
	}
}

#### GREAT tss

for(genome in c("hg19", "hg38", "mm9", "mm10")) {
	gr = randomRegions(genome = genome)

	test_great(gr, "GO:BP", qq("GREAT:@{genome}"))
}


###### gencode

for(v in paste0("v", 18:19)) {
	gr = randomRegions(genome = "hg19")
	test_great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("v", 20:40)) {
	gr = randomRegions(genome = "hg38")
	test_great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("vM", 1)) {
	gr = randomRegions(genome = "mm9")
	test_great(gr, "GO:BP", qq("gencode_@{v}"))
}
for(v in paste0("vM", 2:25)) {
	gr = randomRegions(genome = "mm10")
	test_great(gr, "GO:BP", qq("gencode_@{v}"))
}




##### set background
gr = randomRegions(genome = "hg19")
test_great(gr, "GO:BP", "hg19", background = paste0("chr", 1:22))
test_great(gr, "GO:BP", "hg19", exclude = c("chrX", "chrY"))
test_great(gr, "MSigDB:H", "hg19", exclude = getGapFromUCSC("hg19"))



###### BioMart

wrong_dataset = NULL
for(i in 1:nrow(BIOMART)) {
	dataset = BIOMART[i, 1]

	oe = try({
		gr = randomRegionsFromBioMartGenome(dataset)
		test_great(gr, "GO:BP", biomart_dataset = dataset)
	})
	if(inherits(oe, "try-error")) {
		wrong_dataset = c(wrong_dataset, dataset)
	}
}

wrong_dataset2 = NULL
for(dataset in wrong_dataset) {
	oe = try({
		gr = randomRegionsFromBioMartGenome(dataset)
		test_great(gr, "GO:BP", biomart_dataset = dataset)
	})
	if(inherits(oe, "try-error")) {
		wrong_dataset2 = c(wrong_dataset2, dataset)
	}
}
