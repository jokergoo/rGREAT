context("Test `availableCategories` and `availableOntologies`")

test_that("Test `availableCategories`", {
    
    job = new('GreatJob')
    job@parameters$species = "hg19"
    expect_that(availableCategories(job),
        is_identical_to(c("GO", "Phenotype_data_and_human_desease", "Pathway_Data", "Gene_Expression", "Regulatory_Motifs", "Gene_Families")))
    
    job = new('GreatJob')
    job@parameters$species = "hg18"
    expect_that(availableCategories(job),
        is_identical_to(c("GO", "Phenotype_data_and_human_desease", "Pathway_Data", "Gene_Expression", "Regulatory_Motifs", "Gene_Families")))
    
    job = new('GreatJob')
    job@parameters$species = "mm9"
    expect_that(availableCategories(job),
        is_identical_to(c("GO", "Phenotype_data", "Pathway_Data", "Gene_Expression", "Regulatory_Motifs", "Gene_Families")))

    job = new('GreatJob')
    job@parameters$species = "danRer7"
    expect_that(availableCategories(job),
        is_identical_to(c("GO", "Phenotype_data", "Pathway_Data", "Gene_Expression", "Gene_Families")))

})
