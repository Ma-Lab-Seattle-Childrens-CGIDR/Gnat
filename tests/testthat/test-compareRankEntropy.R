test_that("Compare Rank Entropy Works",{
    testSElist <- generateSummarizedExperiment(
        ngenesOrdered = 10,
        ngenesTotal = 20,
        nsamplesOrdered = 15,
        nsamplesTotal = 30,
        dist = rnorm,
        reorderGenes = FALSE,
        genePrefix = "g_",
        samplePrefix = "s_"
    )
    seObj <- testSElist$seObject
    orderedGenes <- testSElist$orderedGenes
    orderedSamples <- testSElist$orderedSamples
    orderedGeneNames <- testSElist$orderedGeneNames
    orderedSampleNames <- testSElist$orderedSampleNames
    unorderedGenes <- testSElist$unorderedGenes
    unorderedSamples <- testSElist$unorderedSamples
    unorderedGeneNames <- testSElist$unorderedGeneNames
    unorderedSampleNames <- testSElist$unorderedSampleNames
    underlyingMatrix <- SummarizedExperiment::assays(seObj)$geneExpression

    biopparam <- BiocParallel::SerialParam(RNGseed=42)

    expectedResults <- raceComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )

    # Test 1: Test if wrapper works the same as underlying function
    actualResults <- compareRankEntropy(seObj,
                                        method="RACE",
                                        phenotype1=orderedSamples,
                                        phenotype2=unorderedSamples,
                                        geneNetworkList =
                                            list(gn1=orderedGenes,
                                                 gn2=unorderedGenes),
                                        assayName="geneExpression",
                                        bootstrapIterations = 20,
                                        replace=TRUE,
                                        asFrame=TRUE,
                                        BPPARAM = biopparam
                                        )
    expect_equal(actualResults, expectedResults)
    # Test 2: Adjust phenotype and gene input types
    actualResults <- compareRankEntropy(seObj,
                                        method="RACE",
                                        phenotype1=orderedSampleNames,
                                        phenotype2=unorderedSampleNames,
                                        geneNetworkList =
                                            list(gn1=orderedGeneNames,
                                                 gn2=unorderedGeneNames),
                                        assayName="geneExpression",
                                        bootstrapIterations = 20,
                                        replace=TRUE,
                                        asFrame=TRUE,
                                        BPPARAM = biopparam
    )
    expect_equal(actualResults, expectedResults)
    # Test 3: Test Logical Input Types
    actualResults <-
        compareRankEntropy(
            seObj,
            method = "RACE",
            phenotype1 = seObj$phenotypeNum == 1,
            phenotype2 = seObj$phenotypeNum == 2,
            geneNetworkList =
                list(
                    gn1 = SummarizedExperiment::rowData(seObj)$networkNum == 1,
                    gn2 = SummarizedExperiment::rowData(seObj)$networkNum == 2
                ),
            assayName = "geneExpression",
            bootstrapIterations = 20,
            replace = TRUE,
            asFrame = TRUE,
            BPPARAM = biopparam
        )
    expect_equal(actualResults, expectedResults)
    # Test 4: Test if wrapper works with matrix input
    actualResults <- compareRankEntropy(underlyingMatrix,
                                        method="RACE",
                                        phenotype1=orderedSamples,
                                        phenotype2=unorderedSamples,
                                        geneNetworkList =
                                            list(gn1=orderedGenes,
                                                 gn2=unorderedGenes),
                                        assayName="geneExpression",
                                        bootstrapIterations = 20,
                                        replace=TRUE,
                                        asFrame=TRUE,
                                        BPPARAM = biopparam
    )
    expect_equal(actualResults, expectedResults)
    # Test 5: Adjust phenotype and gene input types
    actualResults <- compareRankEntropy(underlyingMatrix,
                                        method="RACE",
                                        phenotype1=orderedSampleNames,
                                        phenotype2=unorderedSampleNames,
                                        geneNetworkList =
                                            list(gn1=orderedGeneNames,
                                                 gn2=unorderedGeneNames),
                                        assayName="geneExpression",
                                        bootstrapIterations = 20,
                                        replace=TRUE,
                                        asFrame=TRUE,
                                        BPPARAM = biopparam
    )
    expect_equal(actualResults, expectedResults)
    # Test 6: Test Logical Input Types
    actualResults <-
        compareRankEntropy(
            underlyingMatrix,
            method = "RACE",
            phenotype1 = seObj$phenotypeNum == 1,
            phenotype2 = seObj$phenotypeNum == 2,
            geneNetworkList =
                list(
                    gn1 = SummarizedExperiment::rowData(seObj)$networkNum == 1,
                    gn2 = SummarizedExperiment::rowData(seObj)$networkNum == 2
                ),
            assayName = "geneExpression",
            bootstrapIterations = 20,
            replace = TRUE,
            asFrame = TRUE,
            BPPARAM = biopparam
        )
    expect_equal(actualResults, expectedResults)
})
