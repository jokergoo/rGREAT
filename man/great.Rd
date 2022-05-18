\name{great}
\alias{great}
\title{
Perform GREAT analysis
}
\description{
Perform GREAT analysis
}
\usage{
great(gr, gene_sets, txdb, chromosomes = NULL,
    mode = "basalPlusExt", basal_upstream = 5000,
    basal_downstream = 1000, extension = 1000000,
    extended_tss = NULL, background = NULL, exclude = NULL,
    verbose = great_opt$verbose)
}
\arguments{

  \item{gr}{A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object. This is the input regions.}
  \item{gene_sets}{A single string of defautly supported gene sets collections (see the full list in "Details" section), or a named list of vectors where each vector correspond to a gene set.}
  \item{txdb}{Name of "TxDb.*" packages from Bioconductor. All supported TxDb packages are in \code{rGREAT:::BIOC_ANNO_PKGS$txdb}. Note short genome version can also be used here such as "hg19" or "hg19.knownGene".}
  \item{chromosomes}{The original TxDb data contains not only the "normal chromosomes" but also a lot of small contigs. This function can automatically identify those "normal chromosomes" (not guaranteed), but users can always control which chromosomes to use by specify a vector of chromosome names to this argument.}
  \item{mode}{The mode to extend genes. Value should be one of 'basalPlusExt', 'twoClosest' and 'oneClosest'. See \code{\link{extendTSS}} for details.}
  \item{basal_upstream}{In 'basalPlusExt' mode, number of base pairs extending to the upstream of TSS to form the basal domains.}
  \item{basal_downstream}{In 'basalPlusExt' mode, number of base pairs extending to the downstream of TSS to form the basal domains.}
  \item{extension}{Extensions from the basal domains.}
  \item{extended_tss}{If your organism is not defaultly supported, you can first prepare one by \code{\link{extendTSSFromDataFrame}} or \code{\link{extendTSS}}, and set the object to this argument.}
  \item{background}{Background regions.}
  \item{exclude}{Regions that are excluded from analysis such as gap regions (which can be get by \code{\link{getGapFromUCSC}}).}
  \item{verbose}{Whether to print messages.}

}
\details{
When \code{background} or \code{exclude} is set, the analysis is restricted in the background regions, still by using Binomial method. Note
this is different from the original GREAT method which uses Fisher's exact test if background regions is set. See \code{\link{submitGreatJob}} for explanations.

rGREAT supports the following GO gene sets for all organisms (note "GO:" can be omitted):

\describe{
  \item{"GO:BP":}{Biological Process, from GO.db package.}
  \item{"GO:CC":}{Cellular Component, from GO.db package.}
  \item{"GO:MP":}{Molecular Function, from GO.db pacakge.}
}

rGREAT also supports gene sets collections from MSigDB (with the msigdbr package, note this is only for human, "msigdb:" can be omitted):

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

If the defaultly supported gene sets and TxDb are used, Entrez gene ID is always used as the main gene ID. If you provide a self-defined
\code{gene_sets} or \code{extended_tss}, you need to make sure they two have the same gene ID types.
}
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
gr = randomRegions()
res = great(gr, "MSigDB:H", "hg19")
}
}
