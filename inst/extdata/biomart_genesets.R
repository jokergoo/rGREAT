

library(GetoptLong)
get_datasets = function(ensembl, label = "") {
	datasets = listDatasets(ensembl)

	write.csv(datasets, file = qq("genesets/00_@{label}.csv"), row.names = FALSE)

	qqcat("####### @{nrow(datasets)} datasets #######\n")
	for(da in datasets[, 1]) {
		qqcat("@{da}\n")

		if(da == "dmelanogaster_eg_gene") next
		
		cat("  gene table...\n")
		at = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id")

		if(!file.exists(qq("genesets/@{da}_genes.rds"))) {
			ensembl = useDataset(dataset = da, mart = ensembl)
			gene = getBM(attributes = at, mart = ensembl)
			saveRDS(gene, file = qq("genesets/@{da}_genes.rds"))
		}

		cat("  go gene sets table...\n")
		at = c("ensembl_gene_id", "go_id", "namespace_1003")

		if(!file.exists(qq("genesets/@{da}_go_genesets.rds"))) {
			ensembl = useDataset(dataset = da, mart = ensembl)
			oe = try({
				go = getBM(attributes = at, mart = ensembl)
				saveRDS(go, file = qq("genesets/@{da}_go_genesets.rds"))
			})

			if(inherits(oe, "try-error")) {
				gene = readRDS(qq("genesets/@{da}_genes.rds"))
				gene_id = gene[, "ensembl_gene_id"]

				go = NULL
				for(i in seq_len(floor(length(gene_id)/1000))) {
					ind = 1:1000 +1000*(i-1)
					ind = ind[ind <= length(gene_id)]
					qqcat("  genes: @{ind[1]} ~ @{ind[length(ind)]}\n")
					gi = gene_id[ind]
					go = rbind(go, getBM(attributes = at, mart = ensembl, filter = "ensembl_gene_id", value = gi))
				}
				saveRDS(go, file = qq("genesets/@{da}_go_genesets.rds"))
			}
		}
	}
}


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", mirror = "uswest")
# attributes = listAttributes(ensembl)
get_datasets(ensembl, "genes_mart")

listEnsemblGenomes()

for(mart in c("protists_mart", "fungi_mart", "metazoa_mart", "plants_mart")) {
	ensembl <- useEnsemblGenomes(biomart = mart)
	get_datasets(ensembl, mart)
}



#### genes saved as GRanges objects:
library(GenomicRanges)
all_files = list.files(path = "genesets", pattern = "genes.rds$", full.names = TRUE)
for(f in all_files) {
	cat(f, "\n")
	df = readRDS(f)
	df = df[df$start_position <= df$end_position, , drop = FALSE]
	gr = GRanges(seqnames = df$chromosome_name, ranges = IRanges(df$start_position, df$end_position), strand = ifelse(df$strand > 0, "+", "-"),
		gene_id = df$ensembl_gene_id)
	saveRDS(gr, qq("genesets_processed/genes/granges_@{basename(f)}"))
}

#### for each GO term, merge all its offspring terms
library(GO.db)

all_files = list.files(path = "genesets", pattern = "go_genesets.rds$", full.names = TRUE)
for(f in all_files) {
	cat(f, "\n")
	df = readRDS(f)
	df = df[df$go_id != "", , drop = FALSE]

	df1 = df[df$namespace_1003 == "biological_process", , drop = FALSE]
	gs = split(df1$ensembl_gene_id, df1$go_id)

	gs2 = lapply(names(gs), function(nm) {
		go_id = c(nm, GOBPOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = names(gs)
	saveRDS(gs2, qq("genesets_processed/genesets/bp_@{basename(f)}"))

	df1 = df[df$namespace_1003 == "cellular_component", , drop = FALSE]
	gs = split(df1$ensembl_gene_id, df1$go_id)

	gs2 = lapply(names(gs), function(nm) {
		go_id = c(nm, GOCCOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = names(gs)
	saveRDS(gs2, qq("genesets_processed/genesets/cc_@{basename(f)}"))

	df1 = df[df$namespace_1003 == "molecular_function", , drop = FALSE]
	gs = split(df1$ensembl_gene_id, df1$go_id)

	gs2 = lapply(names(gs), function(nm) {
		go_id = c(nm, GOMFOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = names(gs)
	saveRDS(gs2, qq("genesets_processed/genesets/mf_@{basename(f)}"))

}


