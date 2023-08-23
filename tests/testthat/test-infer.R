test_that("Entropy calculation works for vectors", {
    testVector <- c(1, 1, 1, 1, 1)
    expect_equal(vectorEntropy(testVector), 0)
    testVector2 <- c(1, 1, 2)
    expect_equal(vectorEntropy(testVector2), 0.9182958,
        tolerance = 1e-4
    )
    testVector3 <- c(1, 2, 3)
    expect_equal(vectorEntropy(testVector3), 1.584963,
        tolerance = 1e-4
    )
})

test_that("Entropy calculation works for matrices", {
    testMatrix <- matrix(c(
        1, 1, 1, 1, 1,
        1, 2, 1, 2, 2,
        1, 2, 3, 4, 5,
        6, 5, 4, 3, 2
    ), ncol = 4)
    expect_length(matrixEntropy(testMatrix, margin = 2), 4)
    expect_length(matrixEntropy(testMatrix, margin = 1), 5)
    expectedEntropy <- c(0, 0.9709506, -log2(1 / 5), -log2(1 / 5))
    expect_equal(
        matrixEntropy(testMatrix, margin = 2),
        expectedEntropy
    )
})

test_that("Infer Rank Function Works", {
    testMat <- matrix(
        c(0.95,2.4,3.2,4.6,6.2,5.7,7.5,
          1.2,1.95,3.95,2.65,5.4,6.2,7.3,
          2.1,1.2,3.3,4.4,5.5,6.4,7.3,
          2.3,3.2,1.1,6.3,4.5,5.3,7.1),
        ncol=4,
        nrow=7
    )
    expRank <- matrix(
        c(1,2,3,4,6,5,7,
          1,2,4,3,5,6,7,
          2,1,3,4,5,6,7,
          2,3,1,6,4,5,7),
        ncol=4,
        nrow=7
    )
    testRank <- inferRankFunction(testMat)
    # Should have same number of rows
    expect_equal(ncol(testRank), ncol(testMat))
    # Should have integer type
    expect_type(testRank, "integer")
    # Should be equal to expRank
    expect_equal(testRank, expRank)
})

test_that("Score Function Works",{
    testMat <- matrix(
        c(0.95,2.4,3.2,4.6,6.2,5.7,7.5,
          1.2,1.95,3.95,2.65,5.4,6.2,7.3,
          2.1,1.2,3.3,4.4,5.5,6.4,7.3,
          2.3,3.2,1.1,6.3,4.5,5.3,7.1),
        ncol=4,
        nrow=7
    )
    rankMat <- inferRankFunction(testMat)
    score <- inferScoreFunction(rankMat)
    expect_type(score, "double")

    identicalRank <- matrix(seq_len(7), ncol=5, nrow=7)
    identicalRankScore <- inferScoreFunction(identicalRank)
    expect_equal(identicalRankScore, 0)

    smallTestMat <- matrix(
        c(1,2,3,
          2,3,1,
          1,2,3,
          1,3,2),
        nrow=3,
        ncol=4
    )
    smallKnownScore <- 1.103759
    expect_equal(inferScoreFunction(smallTestMat), smallKnownScore,
                 tolerance =  1e-4)

})

test_that("Gene Score Works", {

})


test_that("Bootstrap Score Works", {
    # Set seed for cluster
    param = BiocParallel::SnowParam(workers=2, RNGseed = 42)
    # Set seed for datagen
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    testvals <- datagenGenerate(ngenesOrdered = 10,
                                ngenesTotal = 10, nsamplesOrdered = 10,
                                nsamplesTotal = 20, dist = rnorm,
                                reorderGenes=TRUE)
    testExpression <- testvals$expression
    orderedGenes <- testvals$orderedGenes
    orderedSamples <- testvals$orderedSamples
    unorderedSamples <- setdiff(seq_len(20), orderedSamples)
    # Check if known significant results are indeed significant
    actualRes <- inferBootstrapScore(expression=testExpression,
                                     geneNetwork = orderedGenes,
                                     phenotype1 = orderedSamples,
                                     phenotype2 = unorderedSamples,
                                     bootstrapIterations = 100,
                                     replace=TRUE,
                                     BPPARAM=param)
    expect_equal(abs(actualRes$p1Score - actualRes$p2Score),
                 actualRes$absoluteDifference)
    expect_lte(actualRes$pValue, 0.01)
    # Check if randomly ordered phenotypes generate non-significant results
    unordered1 <- sample(unorderedSamples, 5)
    unordered2 <- sample(unorderedSamples, 5)
    randomRes <- inferBootstrapScore(expression=testExpression,
                                     geneNetwork = orderedGenes,
                                     phenotype1 = unordered1,
                                     phenotype2 = unordered2,
                                     bootstrapIterations = 100,
                                     replace=TRUE,
                                     BPPARAM=param)
    expect_gte(randomRes$pValue, 0.1)
    # Test no replacement
    actualRes <- inferBootstrapScore(expression=testExpression,
                                     geneNetwork = orderedGenes,
                                     phenotype1 = orderedSamples,
                                     phenotype2 = unorderedSamples,
                                     bootstrapIterations = 100,
                                     replace=FALSE,
                                     BPPARAM=param)
    expect_no_error(inferBootstrapScore(expression=testExpression,
                                        geneNetwork = orderedGenes,
                                        phenotype1 = orderedSamples,
                                        phenotype2 = unorderedSamples,
                                        bootstrapIterations = 100,
                                        replace=FALSE,
                                        BPPARAM=param))
    expect_lte(actualRes$pValue, 0.05)

})

test_that("Compare Phenotypes Works", {
    # Set seed for cluster
    param <- BiocParallel::SnowParam(workers=2, RNGseed = 42)
    # Set seed for datagen
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    testvals <- datagenGenerate(ngenesOrdered = 15,
                                ngenesTotal = 30, nsamplesOrdered = 10,
                                nsamplesTotal = 20, dist = rnorm,
                                reorderGenes=TRUE)
    testExpression <- testvals$expression
    orderedGenes <- testvals$orderedGenes
    orderedSamples <- testvals$orderedSamples
    unorderedSamples <- setdiff(seq_len(20), orderedSamples)
    # Get subsets of the ordered genes to act as random
    network1 <- sample(orderedGenes, 5)
    network2 <- sample(orderedGenes, 5)
    network3 <- sample(orderedGenes, 5)
    geneNetworkList <- list(network1, network2, network3)
    resActual <- inferComparePhenotypes(expression = testExpression,
                                   geneNetworkList = geneNetworkList,
                                   phenotype1 = orderedSamples,
                                   phenotype2 = unorderedSamples,
                                   bootstrapIterations = 100,
                                   replace=TRUE,asFrame = TRUE,
                                   BPPARAM=param)
    expect_equal(abs(resActual$p1Score-resActual$p2Score),
                 resActual$absoluteDifference)
    expect_true(all(resActual$pValue<0.05))
    # Test unordered samples to see if results are not-significant
    unordered1 <- sample(unorderedSamples, 5)
    unordered2 <- sample(unorderedSamples, 5)
    resUnordered <- inferComparePhenotypes(expression = testExpression,
                                        geneNetworkList = geneNetworkList,
                                        phenotype1 = unordered1,
                                        phenotype2 = unordered2,
                                        bootstrapIterations = 100,
                                        replace=TRUE,asFrame = TRUE,
                                        BPPARAM=param)
    expect_true(all(resUnordered$pValue>0.05))

    # Test list return
    names(geneNetworkList) <- c("A","B","C")
    resList <- inferComparePhenotypes(expression = testExpression,
                                        geneNetworkList = geneNetworkList,
                                        phenotype1 = orderedSamples,
                                        phenotype2 = unorderedSamples,
                                        bootstrapIterations = 100,
                                        replace=TRUE,asFrame = FALSE,
                                        BPPARAM=param)
    expect_equal(names(resList), c("A", "B", "C"))
    for(elem in resList){
        expect_named(elem,
                     c("p1Score", "p2Score", "absoluteDifference", "pValue"))
        expect_equal(abs(elem$p1Score-elem$p2Score), elem$absoluteDifference)
        expect_lte(elem$pValue, 0.05)
    }
})
