# Test Internal Wrapper Functions -----------------------------------------

test_that("Wrap Method Compare Phenotypes Works",{
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
    expectedCraneResults <- craneComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )
    # Test 1
    actualCraneResults <- .wrapMethodComparePhenotypes(
        seObject = seObj,
        method=craneComparePhenotypes,
        phenotype1 = orderedSamples,
        phenotype2 = unorderedSamples,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        assayName = "geneExpression",
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM = biopparam
    )
    expect_equal(actualCraneResults, expectedCraneResults)
    # Test 2
    actualCraneResults <- .wrapMethodComparePhenotypes(
        seObject = seObj,
        method=craneComparePhenotypes,
        phenotype1 = orderedSampleNames,
        phenotype2 = unorderedSampleNames,
        geneNetworkList = list(gn1=orderedGeneNames, gn2=unorderedGeneNames),
        assayName = "geneExpression",
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM = biopparam
    )
    expect_equal(actualCraneResults, expectedCraneResults)
    # Test 3
    actualCraneResults <- .wrapMethodComparePhenotypes(
        seObject = seObj,
        method=craneComparePhenotypes,
        phenotype1 = seObj$phenotypeNum==1,
        phenotype2 = seObj$phenotypeNum==2,
        geneNetworkList = list(gn1=orderedGeneNames, gn2=unorderedGeneNames),
        assayName = "geneExpression",
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM = biopparam
    )
    expect_equal(actualCraneResults, expectedCraneResults)
    # Test 4
    actualCraneResults <- .wrapMethodComparePhenotypes(
        seObject = seObj,
        method=craneComparePhenotypes,
        phenotype1 = seObj$phenotypeNum==1,
        phenotype2 = seObj$phenotypeNum==2,
        geneNetworkList =
            list(gn1=SummarizedExperiment::rowData(seObj)$networkNum==1,
                 gn2=SummarizedExperiment::rowData(seObj)$networkNum==2),
        assayName = "geneExpression",
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM = biopparam
    )
    expect_equal(actualCraneResults, expectedCraneResults)
    # Test 5
    actualCraneResults <- .wrapMethodComparePhenotypes(
        seObject = seObj,
        method=craneComparePhenotypes,
        phenotype1 = seObj$phenotypeStr=="one",
        phenotype2 = seObj$phenotypeStr=="two",
        geneNetworkList =
            list(gn1=SummarizedExperiment::rowData(seObj)$networkStr=="one",
                 gn2=SummarizedExperiment::rowData(seObj)$networkStr=="two"),
        assayName = "geneExpression",
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM = biopparam
    )
    expect_equal(actualCraneResults, expectedCraneResults)

    # Error Test
    expect_error(
        .wrapMethodComparePhenotypes(
            seObject = seObj,
            method=craneComparePhenotypes,
            phenotype1 = seObj$phenotypeStr=="one",
            phenotype2 = seObj$phenotypeStr=="two",
            geneNetworkList =
                list(gn1=SummarizedExperiment::rowData(seObj)$networkStr=="one",
                     gn2=SummarizedExperiment::rowData(seObj)$networkStr=="two"),
            assayName = "counts",
            bootstrapIterations = 20,
            replace=TRUE,
            asFrame=TRUE,
            BPPARAM = biopparam
        ),
        "Assay counts not found in Summarized Experiment"
    )
})

test_that("Wrap Method By Sample Works",{
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

    res <- .wrapMethodBySample(
        seObject = seObj,
        method=craneSampleScore,
        methodName="CRANE",
        phenotype=orderedSamples,
        geneNetwork= orderedGenes,
        colName="sampleEntropy",
        assayName="geneExpression"
    )

    expect_s4_class(res, "SummarizedExperiment")
    expect_true(
        "sampleEntropy" %in% colnames(
            SummarizedExperiment::colData(res)))
    # Since the ordered genes are perfectly ordered, the by sample entropy
    # should all be zero
    sampleEntropyVector <- as.vector(
        SummarizedExperiment::colData(res)[orderedSamples,"sampleEntropy"])
    expect_equal(sampleEntropyVector, rep(0, length(sampleEntropyVector)))

    # Use the unordered genes to calculate some that will not all be 0
    unorderedRes <- .wrapMethodBySample(
        seObject = seObj,
        method=craneSampleScore,
        methodName="CRANE",
        phenotype=orderedSamples,
        geneNetwork= unorderedGenes,
        colName="sampleEntropy",
        assayName="geneExpression"
    )
    sampleEntropyVector <- as.vector(
        SummarizedExperiment::colData(
            unorderedRes)[orderedSamples,"sampleEntropy"])
    names(sampleEntropyVector) <- orderedSampleNames
    # Calculate expected by sample entropy
    sampleEntropyVectorExpected <-  craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples])
    expect_equal(sampleEntropyVector, sampleEntropyVectorExpected)
})


# Test Compare Phenotype Methods ------------------------------------------
test_that("Method Specific Compare Phenotypes Work",{
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
    # CRANE
    expectedCraneResults <- craneComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )
    actualCraneResults <- seCrane(seObject =seObj,
                                  phenotype1=orderedSamples,
                                  phenotype2=unorderedSamples,
                                  geneNetworkList =
                                      list(gn1=orderedGenes,
                                           gn2=unorderedGenes),
                                  assayName = "geneExpression",
                                  bootstrapIterations = 20,
                                  replace = TRUE,
                                  asFrame = TRUE,
                                  BPPARAM=biopparam)
    expect_equal(actualCraneResults, expectedCraneResults)
    # DIRAC
    expectedDiracResults <- diracComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )
    actualDiracResults <- seDirac(seObject =seObj,
                                  phenotype1=orderedSamples,
                                  phenotype2=unorderedSamples,
                                  geneNetworkList =
                                      list(gn1=orderedGenes,
                                           gn2=unorderedGenes),
                                  assayName = "geneExpression",
                                  bootstrapIterations = 20,
                                  replace = TRUE,
                                  asFrame = TRUE,
                                  BPPARAM=biopparam)
    expect_equal(actualDiracResults, expectedDiracResults)
    # RACE
    expectedRaceResults <- raceComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )
    actualRaceResults <- seRace(seObject =seObj,
                                  phenotype1=orderedSamples,
                                  phenotype2=unorderedSamples,
                                  geneNetworkList =
                                      list(gn1=orderedGenes,
                                           gn2=unorderedGenes),
                                  assayName = "geneExpression",
                                  bootstrapIterations = 20,
                                  replace = TRUE,
                                  asFrame = TRUE,
                                  BPPARAM=biopparam)
    expect_equal(actualRaceResults, expectedRaceResults)
    # INFER
    expectedInferResults <- inferComparePhenotypes(
        expression=underlyingMatrix,
        geneNetworkList = list(gn1=orderedGenes, gn2=unorderedGenes),
        phenotype1=orderedSamples,
        phenotype2=unorderedSamples,
        bootstrapIterations = 20,
        replace=TRUE,
        asFrame=TRUE,
        BPPARAM=biopparam
    )
    actualInferResults <- seInfer(seObject =seObj,
                                phenotype1=orderedSamples,
                                phenotype2=unorderedSamples,
                                geneNetworkList =
                                    list(gn1=orderedGenes,
                                         gn2=unorderedGenes),
                                assayName = "geneExpression",
                                bootstrapIterations = 20,
                                replace = TRUE,
                                asFrame = TRUE,
                                BPPARAM=biopparam)
    expect_equal(actualInferResults, expectedInferResults)
})



# Test By Sample Methods ---------------------------------------------

test_that("Method Specific By Sample Methods Work",{
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

    # CRANE
    unorderedRes <- seCraneBySample(seObject=seObj,
                                    phenotype=orderedSamples,
                                    geneNetwork = unorderedGenes,
                                    colName = "craneSampleEntropy",
                                    assayName = "geneExpression"
                                    )
    sampleEntropyVector <- as.vector(
        SummarizedExperiment::colData(
            unorderedRes)[orderedSamples,"craneSampleEntropy"])
    names(sampleEntropyVector) <- orderedSampleNames
    expectedSampleEntropyVector <- craneSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    expect_equal(sampleEntropyVector, expectedSampleEntropyVector)
    # RACE
    unorderedRes <- seRaceBySample(seObject=seObj,
                                    phenotype=orderedSamples,
                                    geneNetwork = unorderedGenes,
                                    colName = "raceSampleEntropy",
                                    assayName = "geneExpression"
    )
    sampleEntropyVector <- as.vector(
        SummarizedExperiment::colData(
            unorderedRes)[orderedSamples,"raceSampleEntropy"])
    names(sampleEntropyVector) <- orderedSampleNames
    expectedSampleEntropyVector <- raceSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    expect_equal(sampleEntropyVector, expectedSampleEntropyVector)
    # DIRAC
    unorderedRes <- seDiracBySample(seObject=seObj,
                                    phenotype=orderedSamples,
                                    geneNetwork = unorderedGenes,
                                    colName = "diracSampleEntropy",
                                    assayName = "geneExpression"
    )
    sampleEntropyVector <- as.vector(
        SummarizedExperiment::colData(
            unorderedRes)[orderedSamples,"diracSampleEntropy"])
    names(sampleEntropyVector) <- orderedSampleNames
    expectedSampleEntropyVector <- diracSampleScore(
        underlyingMatrix[unorderedGenes, orderedSamples]
    )
    expect_equal(sampleEntropyVector, expectedSampleEntropyVector)
})
