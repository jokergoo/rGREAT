all_gmt_files = list.files(path = "~/project/development/rGREAT_genesets/msigdb", pattern = "gmt$", full.names = TRUE)


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

	output = paste0("~/project/development/rGREAT_genesets/msigdb/", basename(f), ".rds")
	saveRDS(lt, file = output, compress = "xz")
}
