test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("Classifier Works",{
    # Create test expression set
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(floor(runif(10 * 8, min = 1, max = 500000)),
                         ncol = 8, nrow = 10)
    geneIndex <- c(1, 4, 5, 7)
    classificationFunction <- .diracEntropyClassifier(expression,
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
