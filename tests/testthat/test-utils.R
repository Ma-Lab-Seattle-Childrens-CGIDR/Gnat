testthat::test_that("Ensure named works", {
    testListUnamed <- list(1, 2, 3, 4)
    testListNamed <- ensureNamed(testListUnamed, prefix = "g_")
    expect_length(testListNamed, length(testListUnamed))
    expect_named(testListNamed, c("g_1", "g_2", "g_3", "g_4"))
    testListNamed <- list(a = 1, b = 2, c = 3, d = 4)
    testListEnsures <- ensureNamed(testListNamed, prefix = "g_")
    expect_length(testListEnsures, length(testListNamed))
    expect_named(testListEnsures, names(testListNamed))
    expect_type(testListEnsures, "list")
})

test_that("Simple Rank Works",{
    testMat <- matrix(runif(20), ncol=4, nrow=5)
    ranked <- simpleRank(testMat)
    expect_equal(nrow(testMat), nrow(ranked))
    expect_equal(ncol(testMat), ncol(ranked))
    orderedMat <- matrix(seq_len(5), ncol=4, nrow=5)
    ranked <- simpleRank(orderedMat)
    expect_equal(ranked, orderedMat)
})

# ADD TESTS FOR ARGUMENT HANDLING FUNCTIONS
