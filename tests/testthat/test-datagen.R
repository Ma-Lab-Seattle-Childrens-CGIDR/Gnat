test_that("Ordered Vector Creation Works", {
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    res <- datagenOrderedVector(16, rnbinom, size = 5, mu = 15)
    # Test that it is the right length
    expect_length(res, 16)
    # Test that it is in order
    expect_equal(res, sort(res))
    # Test another distribution
    resNorm <- datagenOrderedVector(16, rnorm, mean = 10, sd = 2)
    expect_length(resNorm, 16)
    expect_equal(resNorm, sort(resNorm))
    # Test that the results are approximately normal using shapiro-wilkes
    expect_gte(stats::shapiro.test(
        datagenOrderedVector(30, rnorm, mean = 10, sd = 2)
    )$p.value, 0.05)
})

test_that("Unordered vector creation works", {
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    res <- datagenUnorderedVector(16, rnbinom, size = 5, mu = 15)
    expect_length(res, 16)
    expect_true(!all(res == sort(res)))
})

test_that("Ordered matrix creation works", {
    expression <- datagenOrderedMatrix(10, 11, rnbinom, size = 5, mu = 15)
    # Test for correct shape
    expect_equal(nrow(expression), 10)
    expect_equal(ncol(expression), 11)
    # Test that the ordering is true
    resVec <- logical(9)
    for (i in 1:9) {
        resVec[i] <- all(expression[i, ] <= expression[i + 1, ])
    }
    expect_true(all(resVec))
})

test_that("Unordered matrix creation works", {
    expression <- datagenUnorderedMatrix(10, 11, rnbinom, size = 5, mu = 15)
    # Test for correct shape
    expect_equal(nrow(expression), 10)
    expect_equal(ncol(expression), 11)
})

test_that("Generating expression data works", {
    res <- datagenGenerate(
        ngenesOrdered = 10,
        nsamplesOrdered = 12,
        ngenesTotal = 20,
        nsamplesTotal = 18,
        dist = rnbinom, reorderGenes = TRUE,
        size = 5, mu = 15
    )
    expression <- res$expression
    orderedGenes <- res$orderedGenes
    orderedSamples <- res$orderedSamples
    # Test for shape of result
    expect_equal(nrow(expression), 20)
    expect_equal(ncol(expression), 18)
    # Test for length of index vectors
    expect_length(orderedGenes, 10)
    expect_length(orderedSamples, 12)
    # Test that the genes are in fact ordered
    orderedExpression <- expression[orderedGenes, orderedSamples]
    resVec <- logical(9)
    for (i in 1:9) {
        resVec[i] <- all(orderedExpression[i, ] <= orderedExpression[i + 1, ])
    }
    expect_true(all(resVec))
})

test_that("Swap works", {
    expression <- datagenGenerate(
        ngenesOrdered = 10,
        nsamplesOrdered = 12,
        ngenesTotal = 20,
        nsamplesTotal = 18,
        dist = rnbinom, reorderGenes = TRUE,
        size = 5, mu = 15
    )$expression
    res <- datagenSwap(expression = expression, nswaps = 20)
    # Test shape of result
    expect_equal(ncol(res), ncol(expression))
    expect_equal(nrow(res), nrow(expression))
    # Test swapping effects on ordered matrix
    expression <- matrix(1:16, ncol = 4)
    res <- datagenSwap(expression, nswaps = 1)
    # Check that there was exactly one swap
    expect_equal(sum((expression - res) != 0), 2)
})

test_that("Adding noise works", {
    set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
    expression <- matrix(1:36, ncol = 6)
    noisy <- datagenNoise(expression, rnorm, mean = 0, sd = 2)
    # Test shape of result
    expect_equal(ncol(noisy), ncol(expression))
    expect_equal(nrow(noisy), nrow(expression))
    # Test that the noise added is normal
    noise <- as.vector(noisy - expression)
    p.value <- stats::shapiro.test(noise)$p.value
    expect_gte(p.value, 0.05)
})

test_that("Drop out works", {
    expression <- matrix(1:16, ncol = 4)
    res <- datagenDropout(expression, 3 / 16)
    # Test shape of result
    expect_equal(ncol(res), ncol(expression))
    expect_equal(nrow(res), nrow(expression))
    # Test that 3 entries were dropped
    num_drops <- sum((expression - res) != 0)
    expect_equal(num_drops, 3)
})
