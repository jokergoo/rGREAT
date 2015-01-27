\name{availableOntologies}
\alias{availableOntologies}
\title{
  All available ontology names  


}
\description{
  All available ontology names  


}
\usage{
availableOntologies(job, category = NULL)
}
\arguments{
  \item{job}{job object which is returned by \code{\link{submitGreatJob}}}
  \item{category}{categories. All available categories can be get by \code{\link{availableCategories}}}

}
\details{
  All valid values are:  

  \itemize{
    \item GO_Molecular_Function
    \item GO_Biological_Process
    \item GO_Cellular_Component
    \item Mouse_Phenotype
    \item Human_Phenotype
    \item Disease_Ontology
    \item MSigDB_Cancer_Neighborhood
    \item Placenta_Disorders
    \item PANTHER_Pathway
    \item Pathway_Commons
    \item BioCyc_Pathway
    \item MSigDB_Pathway
    \item MGI_Expression_Detected
    \item MSigDB_Perturbation
    \item Transcription_Factor_Targets
    \item MSigDB_Predicted_Promoter_Motifs
    \item MSigDB_miRNA_Motifs
    \item miRNA_Targets
    \item InterPro
    \item TreeFam
    \item HGNC_Gene_Families
    \item Wiki_Pathway
    \item Zebrafish_Wildtype_Expression
    \item Zebrafish_Phenotype
  }

}
\value{
  A vector 


}
\examples{
set.seed(123)
bed = circlize:::generateRandomBed(nr = 1000, nc = 0)
	
#job = submitGreatJob(bed)
load(paste0(system.file(package = "rGREAT"), "/extdata/job.RData"))

availableOntologies(job)
availableOntologies(job, category = "Pathway_Data")
}
