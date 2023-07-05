testthat::test_that("Rank Correlation works for two vectors",{
  x <- c(1,2,3,4,5)
  y <- c(5,4,3,2,1)
  z <- c(1,3,2,4,5)
  expect_equal(rce.rank_correlation.vector(x,x), 1)
  expect_equal(rce.rank_correlation.vector(x,y), -1)
  expect_equal(rce.rank_correlation.vector(x,z), 0.8)
})


testthat::test_that("Rank correlation works for a matrix",{
  test_mat <- matrix(c(1,2,3,4,5,
                       5,4,3,2,1,
                       1,3,2,4,5),
                     byrow=FALSE, ncol=3)
  expected <- c(1,-1,0.8)
  expect_equal(rce.rank_correlation.matrix(test_mat, test_mat[,1]),
               expected)
})

testthat::test_that("Comparing phenotypes using rank correlation entropy works",
                    {
  # Create test random expression matrices
  expression.p1 <- matrix(c(1.9,2.8,3.6,4.6,5.2,
                            2.3,1.1,3.5,4.65,5.7,
                            1.3,3.2,2.1,4.9,5.2,
                            1.1,2.4,3.2,4.6,5.1), ncol=4)
  expression.p2 <- matrix(runif(6*5, min=1, max=5000000), ncol=6)
  expression.combined <- cbind(expression.p1, expression.p2)
  # Test gene index vector
  gene_index <- c(1,3,4)
  # Phenotype index vectors
  p1.index <- c(1,2,3,4)
  p2.index <- c(5,6,7,8,9,10)

  # Warnings are suppressed due to ties, so p-value is approximate
  # Run the rank correlation entropy calculation
  suppressWarnings({result <- rce.compare_phenotypes(
    expression.combined, gene_index, p1.index, p2.index,
    bootstrap_iterations=100, parallel=FALSE, replace=TRUE, seed1=123,
    seed2=456)})
  expect_true(result<0.0001)
  # Repeat, now with parallel operation
  suppressWarnings({result.parallel <- rce.compare_phenotypes(
    expression.combined, gene_index, p1.index, p2.index,
    bootstrap_iterations=100, parallel=TRUE, cores=2, replace=TRUE, seed1=123,
    seed2=456)})
  expect_true(result.parallel<0.0001)

})
