\name{great}
\alias{great}
\title{
Perform GREAT analysis
}
\description{
Perform GREAT analysis
}
\usage{
great(gr, gene_sets, tss_source, biomart_dataset = NULL,
    min_gene_set_size = 5, mode = "basalPlusExt", basal_upstream = 5000,
    basal_downstream = 1000, extension = 1000000,
    extended_tss = NULL, background = NULL, exclude = "gap",
    cores = 1, verbose = great_opt$verbose)
}
\arguments{

  \item{gr}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. This is the input regions. It is important to keep consistent for the chromosome names of the input regions and the internal TSS regions. Use \code{\link{getTSS}} to see the format of internal TSS regions.}
  \item{gene_sets}{A single string of defautly supported gene sets collections (see the full list in "Genesets" section), or a named list of vectors where each vector correspond to a gene set.}
  \item{tss_source}{Source of TSS. See "TSS" section.}
  \item{biomart_dataset}{The value should be in \code{BioMartGOGeneSets::supportedOrganisms}. }
  \item{min_gene_set_size}{Minimal size of gene sets.}
  \item{mode}{The mode to extend genes. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See \code{\link{extendTSS}} for details.}
  \item{basal_upstream}{In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.}
  \item{basal_downstream}{In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.}
  \item{extension}{Extensions from the basal domains.}
  \item{extended_tss}{If your organism is not defaultly supported, you can first prepare one by \code{\link{extendTSSFromDataFrame}} or \code{\link{extendTSS}}, and set the object to this argument. Please see more examples in the vignette.}
  \item{background}{Background regions. The value can also be a vector of chromosome names.}
  \item{exclude}{Regions that are excluded from analysis such as gap regions (which can be get by \code{\link{getGapFromUCSC}}). The value can also be a vector of chromosome names. It also allows a special character value \code{"gap"} so that gap regions for corresponding organism will be removed from the analysis.}
  \item{cores}{Number of cores to use.}
  \item{verbose}{Whether to print messages.}

}
\details{
When \code{background} or \code{exclude} is set, the analysis is restricted in the background regions, still by using Binomial method. Note
this is different from the original GREAT method which uses Fisher's exact test if background regions is set. See \code{\link{submitGreatJob}} for explanations.

By default, gap regions are excluded from the analysis.
}
\section{Tss}{
rGREAT supports TSS from many organisms. The value of \code{tss_source} should be encoded in a special format:

\itemize{
  \item Name of \code{TxDb.*} packages. Supported packages are in \code{rGREAT:::BIOC_ANNO_PKGS$txdb}.
  \item Genome version of the organism, e.g. "hg19". Then the corresponding TxDb will be used.
  \item In a format of \code{RefSeqCurated:$genome} where \code{$genome} is the genome version of an organism, such as hg19. RefSeqCurated subset will be used.
  \item In a format of \code{RefSeqSelect:$genome} where \code{$genome} is the genome version of an organism, such as hg19. RefSeqSelect subset will be used.
  \item In a format of \code{Gencode_v$version} where \code{$version} is gencode version, such as 19 (for human) or M21 for mouse. Gencode protein coding genes will be used.
  \item In a format of \code{GREAT:$genome}, where \code{$genome} can only be mm9, mm10, hg19, hg38. The TSS from GREAT will be used.
}}
\section{Genesets}{
rGREAT supports the following built-in GO gene sets for all organisms (note "GO:" can be omitted):

\describe{
  \item{"GO:BP":}{Biological Process, from GO.db package.}
  \item{"GO:CC":}{Cellular Component, from GO.db package.}
  \item{"GO:MP":}{Molecular Function, from GO.db pacakge.}
}

rGREAT also supports built-in gene sets collections from MSigDB (note this is only for human, "msigdb:" can be omitted):

\describe{
  \item{"msigdb:H"}{Hallmark gene sets.}
  \item{"msigdb:C1"}{Positional gene sets.}
  \item{"msigdb:C2"}{Curated gene sets.}
  \item{"msigdb:C2:CGP"}{C2 subcategory: chemical and genetic perturbations gene sets. }
  \item{"msigdb:C2:CP"}{C2 subcategory: canonical pathways gene sets. }
  \item{"msigdb:C2:CP:BIOCARTA"}{C2 subcategory: BioCarta subset of CP.}
  \item{"msigdb:C2:CP:KEGG"}{C2 subcategory: KEGG subset of CP.}
  \item{"msigdb:C2:CP:PID"}{C2 subcategory: PID subset of CP.}
  \item{"msigdb:C2:CP:REACTOME"}{C2 subcategory: REACTOME subset of CP.}
  \item{"msigdb:C2:CP:WIKIPATHWAYS"}{C2 subcategory: WIKIPATHWAYS subset of CP.}
  \item{"msigdb:C3"}{Regulatory target gene sets.}
  \item{"msigdb:C3:MIR:MIRDB"}{miRDB of microRNA targets gene sets.}
  \item{"msigdb:C3:MIR:MIR_LEGACY"}{MIR_Legacy of MIRDB.}
  \item{"msigdb:C3:TFT:GTRD"}{GTRD transcription factor targets gene sets.}
  \item{"msigdb:C3:TFT:TFT_LEGACY"}{TFT_Legacy.}
  \item{"msigdb:C4"}{Computational gene sets.}
  \item{"msigdb:C4:CGN"}{C4 subcategory: cancer gene neighborhoods gene sets.}
  \item{"msigdb:C4:CM"}{C4 subcategory: cancer modules gene sets.}
  \item{"msigdb:C5"}{Ontology gene sets.}
  \item{"msigdb:C5:GO:BP"}{C5 subcategory: BP subset.}
  \item{"msigdb:C5:GO:CC"}{C5 subcategory: CC subset.}
  \item{"msigdb:C5:GO:MF"}{C5 subcategory: MF subset.}
  \item{"msigdb:C5:HPO"}{C5 subcategory: human phenotype ontology gene sets.}
  \item{"msigdb:C6"}{Oncogenic signature gene sets.}
  \item{"msigdb:C7"}{Immunologic signature gene sets.}
  \item{"msigdb:C7:IMMUNESIGDB"}{ImmuneSigDB subset of C7.}
  \item{"msigdb:C7:VAX"}{C7 subcategory: vaccine response gene sets.}
  \item{"msigdb:C8"}{Cell type signature gene sets.}
}

If the defaultly supported TxDb is used, Entrez gene ID is always used as the main gene ID. If you provide a self-defined
\code{gene_sets} or \code{extended_tss}, you need to make sure they two have the same gene ID types.}
\section{Biomart}{
rGREAT supports a large number of organisms of which the information is retrieved from Ensembl BioMart. The name of a BioMart dataset
can be assigned to argument \code{biomart_dataset}. All supported organisms can be found with \code{BioMartGOGeneSets::supportedOrganisms}.}
\value{
A \code{\link{GreatObject-class}} object. The following methods can be applied on it:

\itemize{
  \item \code{\link{getEnrichmentTable,GreatObject-method}} to retrieve the result table. 
  \item \code{\link{getRegionGeneAssociations,GreatObject-method}} to get the associations between input regions and genes.
  \item \code{\link{plotRegionGeneAssociations,GreatObject-method}} to plot the associations bewteen input regions and genes.
  \item \code{\link{shinyReport,GreatObject-method}} to view the results by a shiny application.
}
}
\examples{
if(FALSE) {
gr = randomRegions(genome = "hg19")
res = great(gr, "MSigDB:H", "txdb:hg19")
res = great(gr, "MSigDB:H", "TxDb.Hsapiens.UCSC.hg19.knownGene")
res = great(gr, "MSigDB:H", "RefSeq:hg19")
res = great(gr, "MSigDB:H", "GREAT:hg19")
res = great(gr, "MSigDB:H", "Gencode_v19")
res = great(gr, "GO:BP", "hsapiens_gene_ensembl")
}
}
