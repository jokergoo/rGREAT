\docType{class}
\name{GREAT_Job-class}
\alias{GREAT_Job}
\alias{GREAT}
\alias{rGREAT}
\alias{GREAT_Job-class}
\alias{getEnrichmentTables}
\alias{availableOntologies}
\alias{availableCategories}
\alias{plotRegionGeneAssociationGraphs}
\title{Class to store and retrieve GREAT results}
\description{
After submitting request to GREAT server, the generated results will be 
available on GREAT server for some time. The \code{GREAT_Job} class is defined 
to store parameters that user has set and result tables what are retrieved 
from GREAT server.
}

\section{Constructor}{
\code{\link{submitGreatJob}} is used to generate a \code{GREAT_Job} instance.
}

\section{Methods}{

Assuming \code{job} is a \code{GREAT_Job} instance which is returned by 
\code{\link{submitGreatJob}}, following methods can be applied:

\describe{
\item{Get enrichment tables from GREAT web server}{
    \preformatted{
getEnrichmentTables(job, ontology = NULL, category = c("GO", "Pathway_Data"),
    request_interval = 30, max_tries = 100)}

    Arguments are:
    \describe{
		\item{\code{job}}{\code{GREAT_Job} instance}
        \item{\code{ontology}}{ontology names. Valid values are in 
            \code{job$availableOntologies()}. \code{ontology} is prior to 
            \code{category} argument.}
        \item{\code{category}}{Pre-defined categories. A category can contain 
            more than one ontologies. Valid values are in 
            \code{job$availableCategories()}}
        \item{\code{request_interval}}{time interval for two requests. Default 
            is 300 seconds.}
        \item{\code{max_tries}}{maximum tries}
    }
        
    Please note there is no FDR column in original tables. Users should 
    calculate by themselves by functions such as \code{\link[stats]{p.adjust}}
        
    The returned value is a list of data frames in which each one corresponds to 
    result for a single ontology. The structure of the data frames are same as 
    the tables available on GREAT website.
}

\item{All available ontology names}{
    \preformatted{
availableOntologies(job, category = NULL)}

    Arguments are:
    \describe{
		\item{\code{job}}{\code{GREAT_Job} instance}
        \item{\code{category}}{one or multiple categories. All available 
            categories can be get by \code{job$availableCategories()}}
    }

    Following ontologies are supported by GREAT: for human (hg19 and hg18)

    "GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component", 
    "Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology", 
    "MSigDB_Cancer_Neighborhood", "Placenta_Disorders", "PANTHER_Pathway", 
    "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway", 
    "MGI_Expression_Detected", "MSigDB_Perturbation", 
    "Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs",
    "MSigDB_miRNA_Motifs", "miRNA_Targets", "InterPro", "TreeFam", 
    "HGNC_Gene_Families".

    For mouse (mm9):

    "GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component", 
    "Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology", 
    "MSigDB_Cancer_Neighborhood", "Placenta_Disorders", "PANTHER_Pathway", 
    "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway", 
    "MGI_Expression_Detected", "MSigDB_Perturbation",
    "Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs",
    "MSigDB_miRNA_Motifs", "miRNA_Targets", "InterPro", "TreeFam".

    For zebrafish (danRer7):

    "GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component", 
    "Wiki_Pathway", "Zebrafish_Wildtype_Expression", "Zebrafish_Phenotype", 
    "InterPro", "TreeFam".

    The returned values is a vector of ontologies.
}


\item{Available categories}{

    \preformatted{
availableCategories(job)}
    
	Arguments are:
    \describe{
		\item{\code{job}}{\code{GREAT_Job} instance}
	}
	
    For human (hg19 and hg18), there are following categories and corresponding
	ontologies:
    
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

    The returned value is a vector of categories.
}


\item{Plot region-gene association figures}{

    \preformatted{
plotRegionGeneAssociationGraphs(job, type = 1:3, ontology = NULL, 
    termID = NULL, request_interval = 30, max_tries = 100)
} 
    Arguments are:
    \describe{
		\item{\code{job}}{\code{GREAT_Job} instance}
        \item{\code{type}}{type of plots, should be in \code{1, 2, 3}}
        \item{\code{ontology}}{ontology name}
        \item{\code{termID}}{term id}
        \item{\code{request_interval}}{time interval for two requests. Default is 300 seconds.}
        \item{\code{max_tries}}{maximum tries}
    }
    
    Generated figures are:  

      \itemize{
        \item association between regions and genes
        \item distribution of distance to TSS
        \item distribution of absolute distance to TSS
      }
      
    If \code{ontology} and \code{termID} are set, only regions and genes corresponding to 
    selected ontology term will be used. Valid value for \code{ontology} is in 
    \code{job$availableOntologies()} and valid value for \code{termID} is from 'id' column 
    in the table which is returned by \code{job$getEnrichmentTables()}.  

     a \code{\link[GenomicRanges]{GRanges}} object. Columns in metadata are:  

    \describe{
      \item{gene}{genes that are associated with corresponding regions}
      \item{distTSS}{distance from the regions to TSS of the associated gene}
    }
      
     The returned values corresponds to whole input regions or only regions in specified ontology term, 
     depending on user's setting. 

}
}}
\author{
  Zuguang gu <z.gu@dkfz.de>  
}
\examples{
\dontrun{
set.seed(123)
bed = circlize::generateRandomBed(nr = 1000, nc =0)
job = submitGreatJob(bed)
tb = getEnrichmentTables(job)
tb = getEnrichmentTables(job, ontology = c("GO_Molecular_Function", "BioCyc_Pathway"))
tb = getEnrichmentTables(job, category = "GO")

availableCategories(job)
availableOntologies(job)
availableOntologies(job, category = "Pathway_Data")

res = plotRegionGeneAssociationGraphs(job)
res = plotRegionGeneAssociationGraphs(job, type = 1)
res = plotRegionGeneAssociationGraphs(job, ontology = "GO_Molecular_Function", 
    termID = "GO:00004984")
}
}
