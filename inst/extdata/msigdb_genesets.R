all_gmt_files = list.files(path = "~/workspace/rGREAT/msigdb", pattern = "gmt$", full.names = TRUE)


read_gmt = function(x) {
	ln = readLines(x)
	lt = strsplit(ln, "\t")

	nm = sapply(lt, function(x) x[1])
	v = lapply(lt, function(x) x[-(1:2)])
	names(v) = nm
	v
}


for(f in all_gmt_files) {
	lt = read_gmt(f)

	output = paste0("~/workspace/rGREAT/genesets_processed/msigdb/", basename(f), ".rds")
	saveRDS(lt, file = output)
}
