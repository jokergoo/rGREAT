[![ Status](https://travis-ci.org/jokergoo/rGREAT.svg)](https://travis-ci.org/jokergoo/rGREAT) [![codecov](https://img.shields.io/codecov/c/github/jokergoo/rGREAT.svg)](https://codecov.io/github/jokergoo/rGREAT) [![bioc](http://www.bioconductor.org/shields/downloads/rGREAT.svg)](http://bioconductor.org/packages/stats/bioc/rGREAT.html) ![bioc](http://www.bioconductor.org/shields/years-in-bioc/rGREAT.svg)

## Client for GREAT Analysis

This package makes [GREAT](http://great.stanford.edu) (Genomic Regions Enrichment of Annotations Tool) 
analysis automatic by constructing a HTTP POST request according to user's input and automatically retrieving results from GREAT web server.

### Install

**rGREAT** is available on Bioconductor (http://bioconductor.org/packages/devel/bioc/html/rGREAT.html)

```r
source("http://bioconductor.org/biocLite.R")
biocLite("rGREAT")
```

If you want the latest version, install it directly from GitHub:

```{r}
library(devtools)
install_github("jokergoo/rGREAT")
```

### Usage

The input data is a `GRanges` object or a _BED_-format data frame, no matter it is sorted or not.

```r
library(rGREAT)
set.seed(123)
bed = circlize::generateRandomBed(nr = 1000, nc = 0)
head(bed)
```

```
##    chr    start      end
## 1 chr1   155726  2608935
## 2 chr1  6134977 10483365
## 3 chr1 11354986 11423447
## 4 chr1 15134641 19115321
## 5 chr1 22692774 23328609
## 6 chr1 23639094 25639077
```

You can get the summary of your job by calling `job` variable.

```r
job = submitGreatJob(bed)
job
```

```
## 
## Submit time: 2014-12-24 20:36:20 
## Species: hg19 
## Background: wholeGenome 
## Model: Basal plus extension 
##   Proximal: 5 kb upstream, 1 kb downstream,
##   plus Distal: up to 1000 kb
## Include curated regulatory domains
## 
## Ontologies enrichment that has been downloaded:
##   None
```

With `job`, we can now retrieve results from GREAT. The first and the primary results are the tables which contain enrichment statistics for the analysis. By default it will retrieve results from three GO Ontologies and all pathway ontologies. All tables contains statistics for all terms no matter they are significant or not. Users can then make filtering by a self-defined cutoff.

```r
tb = getEnrichmentTables(job)
names(tb)
```

```
## [1] "GO_Molecular_Function" "GO_Biological_Process" "GO_Cellular_Component" "PANTHER_Pathway"      
## [5] "Pathway_Commons"       "BioCyc_Pathway"        "MSigDB_Pathway"
```

```r
head(tb[[1]])
```

```
##           ID                                     name Binom_Genome_Fraction Binom_Expected
## 1 GO:0016841                   ammonia-lyase activity          1.117309e-04     0.11228950
## 2 GO:0003796                        lysozyme activity          1.137105e-03     1.14279100
## 3 GO:0004731 purine-nucleoside phosphorylase activity          1.436141e-05     0.01443321
## 4 GO:0008111   alpha-methylacyl-CoA racemase activity          1.837979e-05     0.01847169
## 5 GO:0051635           bacterial cell surface binding          2.087516e-03     2.09795300
## 6 GO:0016034    maleylacetoacetate isomerase activity          2.102749e-05     0.02113263
##   Binom_Observed_Region_Hits Binom_Fold_Enrichment Binom_Region_Set_Coverage Binom_Raw_PValue
## 1                          2             17.811100              0.0019900500      0.005846831
## 2                          5              4.375254              0.0049751240      0.006317406
## 3                          1             69.284640              0.0009950249      0.014329660
## 4                          1             54.136900              0.0009950249      0.018302300
## 5                          6              2.859930              0.0059701490      0.020238010
## 6                          1             47.320190              0.0009950249      0.020911120
##   Hyper_Total_Genes Hyper_Expected Hyper_Observed_Gene_Hits Hyper_Fold_Enrichment
## 1                 5     0.44832060                        2              4.461094
## 2                12     1.07596900                        4              3.717578
## 3                 1     0.08966411                        1             11.152730
## 4                 1     0.08966411                        1             11.152730
## 5                17     1.52429000                        4              2.624173
## 6                 1     0.08966411                        1             11.152730
##   Hyper_Gene_Set_Coverage Hyper_Term_Gene_Coverage Hyper_Raw_PValue
## 1            0.0012570710                0.4000000       0.06690106
## 2            0.0025141420                0.3333333       0.01772786
## 3            0.0006285355                1.0000000       0.08966411
## 4            0.0006285355                1.0000000       0.08966411
## 5            0.0025141420                0.2352941       0.05957762
## 6            0.0006285355                1.0000000       0.08966411
```

Association between genomic regions and genes can be get by `plotRegionGeneAssociationGraphs()`. The function will make the three plots which are same as on GREAT website and returns a `GRanges` object which contains the associations.

```r
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job)
```

![1](https://cloud.githubusercontent.com/assets/449218/5553875/96849e3c-8c33-11e4-8424-b263ee3b6818.png)

```r
res
```

```
## GRanges object with 1716 ranges and 2 metadata columns:
##          seqnames               ranges strand   |        gene   distTSS
##             <Rle>            <IRanges>  <Rle>   | <character> <numeric>
##      [1]     chr1 [  155726,  2608935]      *   |      ATAD3C     -2738
##      [2]     chr1 [ 6134977, 10483365]      *   |      ERRFI1   -222778
##      [3]     chr1 [ 6134977, 10483365]      *   |     SLC45A1    -75219
##      [4]     chr1 [11354986, 11423447]      *   |      PTCHD2   -150078
##      [5]     chr1 [11354986, 11423447]      *   |      UBIAD1     55962
##      ...      ...                  ...    ... ...         ...       ...
##   [1712]     chrY [36594469, 36986560]      *   |        <NA>      <NA>
##   [1713]     chrY [38466369, 40094584]      *   |        <NA>      <NA>
##   [1714]     chrY [44554177, 44654068]      *   |        <NA>      <NA>
##   [1715]     chrY [46572517, 50931029]      *   |        <NA>      <NA>
##   [1716]     chrY [53045400, 56467562]      *   |        <NA>      <NA>
##   -------
##   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

By specifying ontology and term ID, you can get the association in certain term. Here the term ID is from the first column of the data frame which is returned by `getEnrichmentTables()`.

```r
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job, ontology = "GO_Molecular_Function",
    termID = "GO:0004984")
```

![2](https://cloud.githubusercontent.com/assets/449218/5553879/bd7fa216-8c33-11e4-8638-ee3d8348d5c0.png)


### License

GPL (>= 2)
