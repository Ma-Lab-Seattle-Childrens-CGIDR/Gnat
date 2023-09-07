test_that("Dirac Rank Vector Works",{
    testVec1 <- c(1,2,3,4)
    testVec2 <- c(4,3,2,1)
    testVec3 <- c(1,3,2,4)
    expVec1 <- as.logical(c(1,1,1,1,1,1))
    expVec2 <- as.logical(c(0,0,0,0,0,0))
    expVec3 <- as.logical(c(1,1,1,0,1,1))
    expect_equal(diracRankVector(testVec1), expVec1)
    expect_equal(diracRankVector(testVec2), expVec2)
    expect_equal(diracRankVector(testVec3), expVec3)
    expect_error(diracRankVector(c()), "Expression vector must not be empty")
})

test_that("Dirac Rank Function Works",{
    testMat <- matrix(rnorm(50), ncol = 5, nrow=10)
    resMat <- diracRankFunction(testMat)
    expect_equal(ncol(resMat), ncol(testMat))
    expect_equal(nrow(resMat), choose(nrow(testMat), 2))

    testMat2 <- matrix(
        c(1,2,3,4,
          2,1,3,4,
          4,3,2,1,
          1,3,2,4,
          2,2,2,2),
        nrow=4, ncol=5
    )
    expMat <- matrix(
        as.logical(c(1,1,1,1,1,1,
                     0,1,1,1,1,1,
                     0,0,0,0,0,0,
                     1,1,1,0,1,1,
                     0,0,0,0,0,0)
                   ),nrow=6, ncol=5
    )
    resMat2 <- diracRankFunction(testMat2)
    expect_equal(ncol(resMat2), ncol(testMat2))
    expect_equal(nrow(resMat2), choose(nrow(testMat2),2))
    expect_equal(resMat2, expMat)
})

test_that("Dirac Rank Template Works",{
    expression <- matrix(rnorm(30), ncol=6, nrow=5)
    rankMat <- diracRankFunction(expression)
    rankTemplate <- diracRankTemplate(rankMat)
    expect_length(rankTemplate, choose(nrow(expression),2))
    expect_type(rankTemplate, "logical")
    testMat <- matrix(
        c(1,2,4,3,
          1,2,4,3,
          1,3,2,4,
          1,3,2,4,
          1,4,3,2), ncol=5, nrow=4)
    testRank <- diracRankFunction(testMat)
    actRes <- diracRankTemplate(testRank)
    expRes <- as.logical(c(1,1,1,0,1,0))
    expect_equal(actRes, expRes)
})

test_that("Dirac Rank Matching Score Works",{
    testVec1 <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    tempVec1 <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    res1 <- 1.
    expect_equal(diracRankMatchingScore(testVec1,tempVec1), res1)

    testVec2 <- c(FALSE, TRUE, FALSE, TRUE, TRUE)
    tempVec2 <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    res2 <- 0.6
    expect_equal(diracRankMatchingScore(testVec2,tempVec2), res2)

    testVec3 <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    tempVec3 <- c(FALSE, TRUE, FALSE, TRUE, TRUE)
    res3 <- 0.6
    expect_equal(diracRankMatchingScore(testVec3,tempVec3), res3)

    testVec4 <- c(TRUE, TRUE, TRUE, TRUE)
    tempVec4 <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    expect_error(diracRankMatchingScore(testVec4,tempVec4),
                 "rankVector and rankTemplate must be same length")

    expect_error(diracRankMatchingScore(c(), c()),
                 "Neither rankVector, nor rankTemplate can be empty")
})

test_that("Dirac Matrix Matching Score Works",{
    testExp <- matrix(
        c(1,3,2,4,
          4,3,2,1,
          1,3,2,4,
          1,3,2,4,
          1,2,3,4
          ), nrow=4, ncol=5
    )
    testRankMat <- diracRankFunction(testExp)
    testTemplate <- diracRankTemplate(testRankMat)
    matchingScoreActual <- diracMatrixMatchingScore(testRankMat, testTemplate)
    matchingScoreExp <- c(1.,1/6,1.,1.,1.-1/6)
    expect_length(matchingScoreActual, ncol(testExp))
    expect_equal(matchingScoreActual, matchingScoreExp)
})

test_that("Dirac Score Function Works",{
    testExp <- matrix(rnorm(30), ncol=6, nrow=5)
    testRank <- diracRankFunction(testExp)
    testTemplate <- diracRankTemplate(testRank)
    expScore <- mean(diracMatrixMatchingScore(testRank, testTemplate))
    actScore <- diracScoreFunction(testRank)
    expect_equal(actScore, expScore)

    testExp2 <- matrix(1:20, ncol=4, nrow=5)
    testRank2 <- diracRankFunction(testExp2)
    expect_equal(diracScoreFunction(testRank2), 1.)

    testExp3 <- matrix(
        c(1,2,3,4,5,
          5,4,3,2,1), ncol=2, nrow=5
    )
    testRank3 <- diracRankFunction(testExp3)
    expect_equal(diracScoreFunction(testRank3), 0.5)
})

test_that("Dirac Sample Score Works",{
    testExp <- matrix(rnorm(50), nrow=10, ncol=5)
    sampleRes <- diracSampleScore(testExp)
    expect_length(sampleRes, 5)

    testRank <- diracRankFunction(testExp)
    testTemplate <- diracRankTemplate(testRank)
    testMatchingScores <- diracMatrixMatchingScore(rankMatrix = testRank,
                                                   rankTemplate = testTemplate)

    expect_equal(sampleRes, testMatchingScores)

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
    actualRes <- diracBootstrapScore(expression=testExpression,
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
    randomRes <- diracBootstrapScore(expression=testExpression,
                                     geneNetwork = orderedGenes,
                                     phenotype1 = unordered1,
                                     phenotype2 = unordered2,
                                     bootstrapIterations = 100,
                                     replace=TRUE,
                                     BPPARAM=param)
    expect_gte(randomRes$pValue, 0.1)
    # Test no replacement
    actualRes <- diracBootstrapScore(expression=testExpression,
                                     geneNetwork = orderedGenes,
                                     phenotype1 = orderedSamples,
                                     phenotype2 = unorderedSamples,
                                     bootstrapIterations = 100,
                                     replace=FALSE,
                                     BPPARAM=param)
    expect_no_error(diracBootstrapScore(expression=testExpression,
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
    resActual <- diracComparePhenotypes(expression = testExpression,
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
    resUnordered <- diracComparePhenotypes(expression = testExpression,
                                           geneNetworkList = geneNetworkList,
                                           phenotype1 = unordered1,
                                           phenotype2 = unordered2,
                                           bootstrapIterations = 100,
                                           replace=TRUE,asFrame = TRUE,
                                           BPPARAM=param)
    expect_true(all(resUnordered$pValue>0.05))

    # Test list return
    names(geneNetworkList) <- c("A","B","C")
    resList <- diracComparePhenotypes(expression = testExpression,
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

test_that("Classifier Works",{
    # Create test expression set
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)),
                         ncol = 8, nrow = 10)
    geneIndex <- c(1, 4, 5, 7)
    classificationFunction <- diracClassifier(expression,
                                              phenotype1 = c(1,2,3,4),
                                              phenotype2 = c(5,6,7,8),
                                              geneIndex = geneIndex)
    expectedRet <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
    expect_equal(classificationFunction(expression), expectedRet)

    ## Switch where the samples in each phenotype are, and make sure the returns
    ## are still what are expected
    expression <- expression[c(5,1,3,2,6,8,4,7)]
    expectedRet2 <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
})

test_that("Dirac Rank Difference Score Works", {
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(rnorm(8*6), ncol=8, nrow=6)
    p1 <- c(1,2,3,7)
    p2 <- c(4,5,6,8)
    network <- c(1, 2, 3, 5, 6)
    testRes <- diracRankDifferenceScore(expression = expression,
                                        phenotype1 = p1,
                                        phenotype2 = p2,
                                        geneIndex = network)
    expect_length(testRes,2)
    p1RankDiff <- testRes[[1]]
    p2RankDiff <- testRes[[2]]
    expect_length(p1RankDiff, length(p1))
    expect_length(p2RankDiff, length(p2))
    expect_type(p1RankDiff, "double")
    expect_type(p2RankDiff, "double")

    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)), ncol = 8, nrow = 10)
    rankDifferenceScoreExp <- list(
        c(0.5, -0.1666666667, -0.16666667, 0.5),
        c(0 - 0.5, -0.5, -0.166666667, -0.5)
    )
    rankDifferenceScoreAct <- diracRankDifferenceScore(
        expression,
        c(1, 3, 5, 7),
        c(2, 4, 6, 8),
        c(2, 4, 5, 6)
    )
    names(rankDifferenceScoreExp[[1]]) <- c(1,3,5,7)
    names(rankDifferenceScoreExp[[2]]) <- c(2,4,6,8)
    expect_equal(rankDifferenceScoreAct, rankDifferenceScoreExp)
})

test_that("Diract Classification Rate Function",{
    # Create test expression set
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)),
                         ncol = 8, nrow = 10)
    # Create list for gene_index
    gene_index <- c(1, 4, 5, 7, 10)
    # Create expected classificaiton rate variable
    classRateExpected <- 0.875
    # Find the actual classificaiton rate
    classRateActual <- diracClassificationRate(
        gene_index,
        expression,
        c(1, 2, 3, 4),
        c(5, 6, 7, 8)
    )
    expect_type(classRateActual, "double")
    expect_equal(classRateActual, classRateExpected)
})

test_that("Classification rate comparison works", {
    # Create test expression set
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)),
                         ncol = 8, nrow = 10)
    # Create list of gene_index vectors
    geneIndexList <- list(A = c(1, 4, 5), B = c(2, 3, 6), C = c(7, 8, 9, 10))
    # Create expected list
    classRateCompareExp <- c(
        A = 0.750, B = 0.500, C = 0.625
    )
    classRateComparison <- diracClassificationRateCompare(
        expression, c(1, 2, 3, 4), c(5, 6, 7, 8), geneIndexList
    )
    expect_equal(
        classRateComparison,
        classRateCompareExp
    )
})

test_that("Finding Best Classification rate best works", {
    # Create test expression set
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)), ncol = 8,
                         nrow = 10)
    # Create list of gene_index vectors
    geneIndexList <- list(A = c(1, 4, 5), B = c(2, 3, 6), C = c(7, 8, 9, 10))
    bestClassRateExpected <- "A"
    bestClassRateActual <- diracClassificationRateBest(
        expression, c(1, 2, 3, 4), c(5, 6, 7, 8), geneIndexList
    )
    expect_equal(
        bestClassRateActual,
        bestClassRateExpected
    )
})





















