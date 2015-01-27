\docType{class}
\name{GREAT_Job-class}
\alias{GREAT_Job}
\alias{GREAT}
\alias{rGREAT}
\alias{GREAT_Job-class}
\title{Reference class to store and retrieve GREAT results}
\description{
	Reference class to store GREAT results
}

\section{Methods}{

Assuming \code{job} is a \code{GREAT_Job} object, it has following methods:

\describe{
\item{Get enrichment tables from GREAT web server}{
	\preformatted{
job\$getEnrichmentTables(ontology = NULL, category = c("GO", "Pathway_Data"),
    request_interval = 30, max_tries = 100)
}

	Arguments are:
	\describe{
		\item{\code{ontology}}{ontology names. Valid values are in \code{availableOntologies}. \code{ontology} is prior to \code{category} argument.}
		\item{\code{category}}{Pre-defined categories. Valid values are in \code{availableCategories}}
		\item{\code{request_interval}}{interval for two requests. Default is 300 seconds.}
		\item{\code{max_tries}}{maximum tries}
	}
		
	Please note there is no FDR column in original tables. Users should calculate by themselves
	by functions such as \code{\link[stats]{p.adjust}}
		
	a list of data frames which are as same as in GREAT website.
}

\item{All available ontology names}{
	\preformatted{
job\$availableOntologies(category = NULL)
	}

	Arguments are:
	\describe{
		\item{\code{category}}{categories. All available categories can be get by \code{availableCategories}}
	}
	
	All valid values are: "GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component", "Mouse_Phenotype",
	"Human_Phenotype", "Disease_Ontology", "MSigDB_Cancer_Neighborhood", "Placenta_Disorders",
	"PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway", "MGI_Expression_Detected",
	"MSigDB_Perturbation", "Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs",
	"MSigDB_miRNA_Motifs", "miRNA_Targets", "InterPro", "TreeFam", "HGNC_Gene_Families",
	"Wiki_Pathway", "Zebrafish_Wildtype_Expression", "Zebrafish_Phenotype"

}


\item{available categories}{

	\preformatted{
job\$availableCategories()
}
	
	For human (hg19 and hg18):
	
	\describe{
		\item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
		\item{Phenotype_data_and_human_desease}{"Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology", "MSigDB_Cancer_Neighborhood", "Placenta_Disorders"}
		\item{Pathway_Data}{"PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"}
		\item{Gene_Expression}{"MGI_Expression_Detected", "MSigDB_Perturbation"}
		\item{Regulatory_Motifs}{"Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"}
		\item{Gene_Families}{"InterPro", "TreeFam", "HGNC_Gene_Families"}
	}
	
	For mouse (mm9):
	
	\describe{
		\item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
		\item{Phenotype_data}{"Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology"}
		\item{Pathway_Data}{"PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"}
		\item{Gene_Expression}{"MGI_Expression_Detected", "MSigDB_Perturbation"}
		\item{Regulatory_Motifs}{"Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"}
		\item{Gene_Families}{"InterPro", "TreeFam"}
	}
	
	For zebrafish (danRer7):
	
	\describe{
		\item{GO}{"GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"}
		\item{Phenotype_data}{"Zebrafish_Phenotype"}
		\item{Pathway_Data}{"Wiki_Pathway"}
		\item{Gene_Expression}{"Zebrafish_Wildtype_Expression"}
		\item{Gene_Families}{"InterPro", "TreeFam"}
	}
}


\item{Plot region-gene association figures}{

	\preformatted{
job\$plotRegionGeneAssociationGraphs(type = 1:3, ontology = NULL, 
    termID = NULL, request_interval = 30, max_tries = 100)
} 
	Arguments are:
	\describe{
		\item{\code{type}}{type of plots, should be in \code{1, 2, 3}}
		\item{\code{ontology}}{ontology name}
		\item{\code{termID}}{term id}
		\item{\code{request_interval}}{interval for two requests. Default is 300 seconds.}
		\item{\code{max_tries}}{maximum tries}
	}
	
	Figures are:  

	  \itemize{
		\item association between regions and genes
		\item distribution of distance to TSS
		\item distribution of absolute distance to TSS
	  }
	  
	If \code{ontology} and \code{termID} are set, only regions and genes corresponding to 
	selected ontology term will be used. Valid value for \code{ontology} is in 
	\code{availableOntologies} and valid value for \code{termID} is from 'id' column 
	in the table which is returned by \code{getGreatTable}.  

	 a \code{GRanges} object. Columns in metadata are:  

	\describe{
	  \item{gene}{genes that are associated with corresponding regions}
	  \item{distTSS}{distance from the regions to TSS of the associated gene}
	}
	  
	  The returned values corresponds to whole input regions or only regions in specified ontology term, depending on user's settings. 

}
}}
\author{
  Zuguang gu <z.gu@dkfz.de>  
}
\examples{

}
