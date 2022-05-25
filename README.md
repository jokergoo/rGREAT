# GREAT Analysis - Functional Enrichment on Genomic Regions

[![R-CMD-check](https://github.com/jokergoo/rGREAT/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/rGREAT/actions)
[![codecov](https://img.shields.io/codecov/c/github/jokergoo/rGREAT.svg)](https://codecov.io/github/jokergoo/rGREAT)
[![bioc](https://bioconductor.org/shields/downloads/devel/rGREAT.svg)](https://bioconductor.org/packages/stats/bioc/rGREAT/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/rGREAT.svg)](http://bioconductor.org/packages/devel/bioc/html/rGREAT.html)


**GREAT** ([Genomic Regions Enrichment of Annotations Tool](http://great.stanford.edu)) is a type of
functional enrichment analysis directly performed on genomic regions. This package 
implements the GREAT algorithm (the local GREAT analysis), also it supports directly 
interacting with the GREAT web service (the online GREAT analysis). Both analysis 
can be viewed by a Shiny application.

## Install

**rGREAT** is available on Bioconductor (http://bioconductor.org/packages/devel/bioc/html/rGREAT.html)

```r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("rGREAT")
```

If you want the latest version, install it directly from GitHub:

```r
library(devtools)
install_github("jokergoo/rGREAT")
```

## Online GREAT analysis

With online GREAT analysis, the input regions will be directly submitted to GREAT server, and the results
are automatically retrieved from GREAT server.

```r
set.seed(123)
gr = randomRegions(nr = 1000, genome = "hg19")

job = submitGreatJob(gr)
tbl = getEnrichmentTables(job)
```

## Local GREAT analysis

**rGREAT** also implements the GREAT algorithms locally and it can be seamlessly integrated
to the Bioconductor annotation ecosystem. This means, theoretically, with **rGREAT**, it is possible to perform GREAT analysis
with any organism and with any type of gene set collection / ontology

```r
res = great(gr, "MSigDB:H", "TxDb.Hsapiens.UCSC.hg19.knownGene")
tb = getEnrichmentTable(res)
```

To apply `great()` on other organisms, set the `biomart_dataset` argument:

```r
# giant panda
great(gr, "GO:BP", biomart_dataset = "amelanoleuca_gene_ensembl")
```

For more details, please go to [the package vignettes](https://jokergoo.github.io/rGREAT/).

## License

MIT @ Zuguang Gu
