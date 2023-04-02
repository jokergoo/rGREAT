library(rGREAT)
set.seed(123)
gr = randomRegions(nr = 1000, genome = "hg19")
job = submitGreatJob(gr)
tbl = getEnrichmentTables(job)
plotRegionGeneAssociations(job)
getRegionGeneAssociations(job)
plotRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")
getRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")

saveRDS(job, file = "~/project/development/rGREAT/inst/extdata/GreatJob.rds", compress = "xz")
