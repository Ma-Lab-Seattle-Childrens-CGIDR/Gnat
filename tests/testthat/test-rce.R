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

testthat::test_that("Computing rce for one sample works",{
  # Create test expression matrix
  expression <- matrix(c(1.9,2.8,3.6,4.6,5.2,
                            2.3,1.1,3.5,4.65,5.7,
                            1.3,3.2,2.1,4.9,5.2,
                            1.1,2.4,3.2,4.6,5.1), ncol=4)
  mean_rce.actual <- rce.sample(3, expression)
  mean_rce.expected <- 0.8666667
  expect_true(is.double(mean_rce.actual))
  expect_equal(mean_rce.actual, mean_rce.actual, tolerance=1e-05)
})


testthat::test_that("Creating rce distribution works",{
  # Create test random expression matrices
  expression.p1 <- matrix(c(1.9,2.8,3.6,4.6,5.2,
                            2.3,1.1,3.5,4.65,5.7,
                            1.3,3.2,2.1,4.9,5.2,
                            1.1,2.4,3.2,4.6,5.1), ncol=4)
  expression.p2 <- matrix(runif(6*5, min=1, max=5000000), ncol=6)
  expression.combined <- cbind(expression.p1, expression.p2)
  # gene index vector
  gene_index <- c(1,2,3,5)
  # phenotype vector
  phenotype <- c(1,2,3,4)
  # number of iterations to perform
  iterations=20
  # Test serial operation
  dist <- rce.distribution.central(gene_expression = expression.combined,
                                   gene_index = gene_index,
                                   phenotype=phenotype,
                                   iterations=iterations,
                                   parallel = FALSE,
                                   replace=TRUE,
                                   seed=123)
  # Make sure the distribution is the right length
  expect_vector(dist, size=iterations)
  expect_type(dist, "double")
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
    gene_index, expression.combined, p1.index, p2.index,
    bootstrap_iterations=100, parallel=FALSE, replace=TRUE, seed1=123,
    seed2=456)})
  expect_true(result$p.value<0.0001)
  # Repeat, now with parallel operation
  suppressWarnings({result.parallel <- rce.compare_phenotypes(
    gene_index, expression.combined, p1.index, p2.index,
    bootstrap_iterations=100, parallel=TRUE, cores=2, replace=TRUE, seed1=123,
    seed2=456)})
  expect_true(result.parallel$p.value<0.0001)
})
