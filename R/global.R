
BASE_URL = "http://bejerano.stanford.edu/great/public/cgi-bin"


# the values will be used as names as a list, so underline is used instead of blank space
ONTOLOGY_KEYS = c(
    "GO_Molecular_Function" = "GOMolecularFunction",
    "GO_Biological_Process" = "GOBiologicalProcess",
    "GO_Cellular_Component" = "GOCellularComponent",
    "Mouse_Phenotype" = "MGIPhenotype",
    "Human_Phenotype" = "HumanPhenotypeOntology",
    "Disease_Ontology" = "OsborneAnnotatedDiseaseOntology",
    "MSigDB_Cancer_Neighborhood" = "MSigDBGeneSetsCancerNeighborhood",
    "Placenta_Disorders" = "PlacentaDisorders",
    "PANTHER_Pathway" = "panther",
    "Pathway_Commons" = "pathway",
    "BioCyc_Pathway" = "BioCycPathway",
    "MSigDB_Pathway" = "MSigDBGeneSetsCanonicalPathway",
    "MGI_Expression_Detected" = "MGIExpressionDetected",
    "MSigDB_Perturbation" = "MSigDBGeneSetsPerturbation",
    "Transcription_Factor_Targets" = "TFregulatedGenes",
    "MSigDB_Predicted_Promoter_Motifs" = "MSigDBGeneSetsPromoterMotifs",
    "MSigDB_miRNA_Motifs" = "MSigDBGeneSetsMicroRNAMotifs",
    "miRNA_Targets" = "miRNAregulatedGenes",
    "InterPro" = "interpro",
    "TreeFam" = "treefam",
    "HGNC_Gene_Families" = "geneFamilies",

    # last three are for zebrafish only
    "Wiki_Pathway" = "wikiPathway",
    "Zebrafish_Wildtype_Expression" = "expressionWt",
    "Zebrafish_Phenotype" = "phenotype"
)

ONTOLOGY_NAMES = names(ONTOLOGY_KEYS)

URL_TEMPLATE = NULL
for(nm in ONTOLOGY_NAMES) {
    URL_TEMPLATE[nm] = qq("`BASE_URL`/readJsFromFile.php?path=/scratch/great/tmp/results/@{jobid}.d/`ONTOLOGY_KEYS[nm]`.js", code.pattern = "`CODE`")
}

## pre-defined category for ontologies
CATEGORY = list(
    hg19 = list(
        "GO"                               = c("GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"),
        "Phenotype_data_and_human_desease" = c("Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology", "MSigDB_Cancer_Neighborhood", "Placenta_Disorders"),
        "Pathway_Data"                     = c("PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"),
        "Gene_Expression"                  = c("MGI_Expression_Detected", "MSigDB_Perturbation"),
        "Regulatory_Motifs"                = c("Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"),
        "Gene_Families"                    = c("InterPro", "TreeFam", "HGNC_Gene_Families")
    ),
    mm9 = list(
        "GO"                               = c("GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"),
        "Phenotype_data"                   = c("Mouse_Phenotype", "Human_Phenotype", "Disease_Ontology"),
        "Pathway_Data"                     = c("PANTHER_Pathway", "Pathway_Commons", "BioCyc_Pathway", "MSigDB_Pathway"),
        "Gene_Expression"                  = c("MGI_Expression_Detected", "MSigDB_Perturbation"),
        "Regulatory_Motifs"                = c("Transcription_Factor_Targets", "MSigDB_Predicted_Promoter_Motifs", "MSigDB_miRNA_Motifs", "miRNA_Targets"),
        "Gene_Families"                    = c("InterPro", "TreeFam")
    ),
    danRer7 = list(
        "GO"                               = c("GO_Molecular_Function", "GO_Biological_Process", "GO_Cellular_Component"),
        "Phenotype_data"                   = c("Zebrafish_Phenotype"),
        "Pathway_Data"                     = c("Wiki_Pathway"),
        "Gene_Expression"                  = c("Zebrafish_Wildtype_Expression"),
        "Gene_Families"                    = c("InterPro", "TreeFam")
    )
)
CATEGORY$hg18 = CATEGORY$hg19

rGREAT_env = new.env()
rGREAT_env$LAST_REQUEST_TIME = 0

