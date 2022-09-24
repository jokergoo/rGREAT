
if(0) {

	
GREAT.options(sleep.time = 60)

set.seed(123)
gr = circlize:::generateRandomBed(nr = 1000, nc = 0)
gr2 = gr
gr2[[1]] = gsub("chr", "", gr2[[1]])
job = submitGreatJob(gr2)
job = submitGreatJob(gr)

job

availableCategories(job)
availableOntologies(job)
availableOntologies(job, category = "GO")
availableOntologies(job, category = c("GO", "Pathway_Data"))

tb = getGreatTable(job)
tb = getGreatTable(job, ontology = c("GO_Molecular_Function", "BioCyc_Pathway"))
tb = getGreatTable(job, category = c("GO", "Pathway_Data"))

par(mfrow = c(3, 1))
df = plotRegionGeneAssociationGraphs(job)

df = plotRegionGeneAssociationGraphs(job, ontology = "GO_Molecular_Function", termID = "GO:0004984")

# test other species
gr = circlize:::generateRandomBed(nr = 1000, nc = 0, species = "hg18")
job = submitGreatJob(gr)
job = submitGreatJob(gr, species = "hg18")
availableCategories(job)
availableOntologies(job)
tb = getGreatTable(job)
df = plotRegionGeneAssociationGraphs(job)
df = plotRegionGeneAssociationGraphs(job, ontology = "GO_Molecular_Function", termID = "GO:0004984")

gr = circlize:::generateRandomBed(nr = 1000, nc = 0, species = "mm9")
job = submitGreatJob(gr)
job = submitGreatJob(gr, species = "mm9")
availableCategories(job)
availableOntologies(job)
tb = getGreatTable(job)
df = plotRegionGeneAssociationGraphs(job)
df = plotRegionGeneAssociationGraphs(job, ontology = "GO_Molecular_Function", termID = "GO:0004984")


}

genes = getGenesFromGencode("v19")
genes = genes[seqnames(genes) == "chr1"]
strand(genes) = "*"
gr = sort(genes)

test_that("test reduce and reduce_by_start_and_end", {
	expect_that(rGREAT:::reduce_by_start_and_end(start(gr), end(gr)), equals(sum(width(reduce(gr)))))
})
