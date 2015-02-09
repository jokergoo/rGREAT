\name{availableCategories-GreatJob-method}
\alias{availableCategories}
\alias{availableCategories,GreatJob-method}
\title{
Available ontology categories  


}
\description{
Available ontology categories  


}
\usage{
\S4method{availableCategories}{GreatJob}(job)
}
\arguments{

  \item{job}{a \code{\link{GreatJob}} instance}

}
\details{
On GREAT, There are several pre-defined ontology categories:  

For human (hg19 and hg18), there are following categories and corresponding ontologies:  

\describe{
  \item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
  \item{Phenotype_data_and_human_desease}{"Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology", "MSigDB_Cancer_Neighborhood", "Placenta_Disorders"}
  \item{Pathway_Data}{"PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"}
  \item{Gene_Expression}{"MGI_Expression_Detected", "MSigDB_Perturbation"}
  \item{Regulatory_Motifs}{"Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"}
  \item{Gene_Families}{"InterPro", "TreeFam", "HGNC_Gene_Families"}
}

For mouse (mm9):  

\describe{
  \item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
  \item{Phenotype_data}{"Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology"}
  \item{Pathway_Data}{"PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"}
  \item{Gene_Expression}{"MGI_Expression_Detected", "MSigDB_Perturbation"}
  \item{Regulatory_Motifs}{"Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"}
  \item{Gene_Families}{"InterPro", "TreeFam"}
}

For zebrafish (danRer7):  

\describe{
  \item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
  \item{Phenotype_data}{"Zebrafish_Phenotype"}
  \item{Pathway_Data}{"Wiki_Pathway"}
  \item{Gene_Expression}{"Zebrafish_Wildtype_Expression"}
  \item{Gene_Families}{"InterPro", "TreeFam"}
}


}
\value{
The returned value is a vector of categories.  


}
\author{
Zuguang gu <z.gu@dkfz.de>  


}
\examples{
job = readRDS(paste0(system.file("extdata", package = "rGREAT"), "/job.rds"))
availableCategories(job)
}
