test_that("By Sample Entropy Works",{
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
    # Test 1: Test that wrapper works as expected
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    seResults <- bySampleEntropy(
        expression=seObj,
        method="CRANE",
        phenotype=orderedSamples,
        geneNetwork = unorderedGenes,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    actualSampleEntropy <- as.vector(
        SummarizedExperiment::colData(
            seResults)[orderedSampleNames,"sampleEntropy"]
        )
    names(actualSampleEntropy) <- orderedSampleNames
    expect_equal(actualSampleEntropy, expectedResults)
    # Test 2: Adjust phenotype/gene network input types
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    seResults <- bySampleEntropy(
        expression=seObj,
        method="CRANE",
        phenotype=orderedSampleNames,
        geneNetwork = unorderedGeneNames,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    actualSampleEntropy <- as.vector(
        SummarizedExperiment::colData(
            seResults)[orderedSampleNames,"sampleEntropy"]
    )
    names(actualSampleEntropy) <- orderedSampleNames
    expect_equal(actualSampleEntropy, expectedResults)
    # Test 3: Adjust phenotype/gene network input types
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    seResults <- bySampleEntropy(
        expression=seObj,
        method="CRANE",
        phenotype=seObj$phenotypeNum==1,
        geneNetwork = SummarizedExperiment::rowData(seObj)$networkNum==2,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    actualSampleEntropy <- as.vector(
        SummarizedExperiment::colData(
            seResults)[orderedSampleNames,"sampleEntropy"]
    )
    names(actualSampleEntropy) <- orderedSampleNames
    expect_equal(actualSampleEntropy, expectedResults)
    # Test 4: Test array input
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    actualResults <- bySampleEntropy(
        expression=underlyingMatrix,
        method="CRANE",
        phenotype=orderedSamples,
        geneNetwork = unorderedGenes,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    expect_equal(actualResults, expectedResults)
    # Test 5: Adjust gene network/phenotype input type
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    actualResults <- bySampleEntropy(
        expression=underlyingMatrix,
        method="CRANE",
        phenotype=orderedSampleNames,
        geneNetwork = unorderedGeneNames,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    expect_equal(actualResults, expectedResults)
    # Test 6: Try Logical vector input
    expectedResults <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    actualResults <- bySampleEntropy(
        expression=underlyingMatrix,
        method="CRANE",
        phenotype=seObj$phenotypeNum==1,
        geneNetwork = SummarizedExperiment::rowData(seObj)$networkNum==2,
        colDataName="sampleEntropy",
        assayName = "geneExpression"
    )
    expect_equal(actualResults, expectedResults)
})
