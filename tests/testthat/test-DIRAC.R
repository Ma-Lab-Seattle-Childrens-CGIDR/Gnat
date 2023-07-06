# Main Function Tests -----------------------------------------------------
testthat::test_that("Compare Network Classification Works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create list of gene_index vectors
  gene_index_list <- list(A=c(1,4,5), B=c(2,3,6), C=c(7,8,9,10))
  # Create expected list
  classification_rate_comparison.expected <- dirac.classification_rate.compare(
    expression, c(1,2,3,4), c(5,6,7,8), gene_index_list, parallel=FALSE,
    as.frame=FALSE)
  classification_rate_comparison.serial <-
    DIRAC.compare_network_classification(expression, c(1,2,3,4), c(5,6,7,8),
                                         gene_index_list, parallel=FALSE,
                                         as.frame=FALSE)
  classification_rate_comparison.parallel <-
    DIRAC.compare_network_classification(expression, c(1,2,3,4), c(5,6,7,8),
                                         gene_index_list, parallel=TRUE,
                                         cores=2,as.frame=FALSE)
  expect_equal(classification_rate_comparison.serial,
               classification_rate_comparison.expected)
  expect_equal(classification_rate_comparison.parallel,
               classification_rate_comparison.expected)
  classification_rate_comparison.frame.expected <- data.frame(
    gene_network=c("A","B","C"), classification_rate=c(0.75, 0.5, 0.625)
  )
  rownames(classification_rate_comparison.frame.expected) <- c("A","B","C")
  classification_rate_comparison.frame.actual <-
    DIRAC.compare_network_classification(expression, c(1,2,3,4), c(5,6,7,8),
                                         gene_index_list, parallel=TRUE,
                                         cores=2, as.frame=TRUE)
  expect_equal(classification_rate_comparison.frame.actual,
               classification_rate_comparison.frame.expected)

})

testthat::test_that("Compare phenotypes works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create list of gene_index vectors
  gene_index_list <- list(A=c(1,4,5), B=c(2,3,6), C=c(7,8,9,10))
  # Test parallel operation
  phenotype_comp.expected <- data.frame(
    gene_network=c("A","B","C"),
    value=c(0,0,0),
    pval=c(0.85,0.65,0.9)
  )
  rownames(phenotype_comp.expected) <- c("A","B","C")
  phenotype_comp.actual <-
    DIRAC.compare_phenotypes(expression, c(1,2,3,4), c(5,6,7,8),
                             gene_index_list, bootstrap_iterations = 20,
                             parallel=TRUE, cores=2, replace=TRUE, seed = 42,
                             as.frame=TRUE)
  expect_equal(phenotype_comp.actual, phenotype_comp.expected)
  # Test serial operation
  phenotype_comp.serial.expected <- data.frame(
    gene_network=c("A","B","C"),
    value=c(0,0,0),
    pval=c(0.7,0.8,0.85)
  )
  rownames(phenotype_comp.serial.expected) <- c("A","B","C")
  phenotype_comp.serial.actual <-
    DIRAC.compare_phenotypes(expression, c(1,2,3,4), c(5,6,7,8),
                             gene_index_list, bootstrap_iterations = 20,
                             parallel=FALSE, cores=2, replace=TRUE, seed = 42,
                             as.frame=TRUE)
  expect_equal(phenotype_comp.serial.actual, phenotype_comp.serial.expected)
})

testthat::test_that("Creating a classifier works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create gene_index vector
  gene_index <- c(1,4,5,7)
  # Create classification functions
  class_func1 <- DIRAC.classifier(expression, c(1,2,3,4), c(5,6,7,8),
                                  gene_index = gene_index,
                                  return="is phenotype1",
                                  names=c("phenotype1", "phenotype2"))
  class_func2 <- DIRAC.classifier(expression, c(1,2,3,4), c(5,6,7,8),
                                  gene_index = gene_index,
                                  return="is phenotype2",
                                  names=c("phenotype1", "phenotype2"))
  class_func3 <- DIRAC.classifier(expression, c(1,2,3,4), c(5,6,7,8),
                                  gene_index = gene_index,
                                  return="name",
                                  names=c("A", "B"))
  class_func4 <- DIRAC.classifier(expression, c(1,2,3,4), c(5,6,7,8),
                                  gene_index = gene_index,
                                  return="phenotype number",
                                  names=c("phenotype1", "phenotype2"))
  class_func5 <- DIRAC.classifier(expression, c(1,2,3,4), c(5,6,7,8),
                                  gene_index = gene_index,
                                  return="rank matching diff",
                                  names=c("phenotype1", "phenotype2"))
  # Vectors with expected return values
  expected.return1 <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
  expected.return2 <- !expected.return1
  expected.return3 <- c("A", "A", "A", "A", "B", "A", "B", "B")
  expected.return4 <- c(1,1,1,1,2,1,2,2)
  expected.return5 <- c(0.5000000,  0.1666667,  0.1666667,  0.5000000,
                        -0.1666667,  0.1666667, -0.1666667, -0.5000000)
  # Vectors with actual return values
  actual.return1 <- class_func1(expression)
  actual.return2 <- class_func2(expression)
  actual.return3 <- class_func3(expression)
  actual.return4 <- class_func4(expression)
  actual.return5 <- class_func5(expression)
  # Test that the functions return the correct values
  expect_equal(actual.return1, expected.return1)
  expect_equal(actual.return2, expected.return2)
  expect_equal(actual.return3, expected.return3)
  expect_equal(actual.return4, expected.return4)
  expect_equal(actual.return5, expected.return5, tolerance=1e-05)
})

# Helper Function Tests ---------------------------------------------------
testthat::test_that("Rank vector Works",{
  # Create expression vector for testing
  expression=c(4,2,1,3)
  # Test if creates correct rank vector for these expression values
  expect_equal(dirac.rank_vector(expression),
               c(FALSE,FALSE,FALSE,FALSE, TRUE, TRUE))
})

testthat::test_that("Rank matrix works",{
  # Create expression matrix for testing
  expression=matrix(c(1,4,3,2,5,4,2,3,1,6,6,1,5,7,9,16,10,5,4,3),
                    nrow = 4, ncol = 5, byrow = TRUE)
  # Run the dirac.rank_matrix function
  rank_matrix <- dirac.rank_matrix(expression = expression)
  # Test that all the rank vectors are correct
  for(col in 1:ncol(expression)){
    expression_vector <- expression[,col]
    rank_vector.expected <- dirac.rank_vector(expression_vector)
    rank_vector.actual <- rank_matrix[,col]
    expect_equal(rank_vector.actual, rank_vector.expected)
  }
})

testthat::test_that("Rank Template Works",{
  # Create expression matrix
  expression=matrix(c(1,4,3,2,5,4,2,3,1,6,6,1,5,7,9,16,10,5,4,3),
                    nrow = 4, ncol = 5, byrow = TRUE)
  template.expected <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
  template.actual <- dirac.rank_template(dirac.rank_matrix(expression))
  expect_equal(template.actual, template.expected)
})

testthat::test_that("Rank matching score works",{
  rank_vector <- c(1,0,1,1,0)
  rank_template <- c(0,1,1,1,0)
  rank_matching_score <- dirac.rank_matching_score(
    rank_vector = rank_vector,
    rank_template = rank_template)
  expect_equal(rank_matching_score, 0.6)
})


testthat::test_that("Rank matching score for a matrix works",{
  rank_template <- c(1,0,1,1,0)
  rank_matrix <- matrix(c(1,0,1,1,0,
                          0,0,1,1,0,
                          1,1,1,1,1,
                          0,0,1,0,0), nrow=5)
  rank_matching_score.expected <- c(1.0,0.8,0.6,0.6)
  rank_matching_score.test <- dirac.rank_matching_score.vector(rank_matrix,
                                                               rank_template)
  expect_equal(rank_matching_score.test, rank_matching_score.expected)
})

testthat::test_that("Rank difference score works",{
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  rank_difference_score.expected <- list(
    c(0.5, -0.1666666667, -0.16666667, 0.5),
    c(0-0.5,-0.5,-0.166666667, -0.5))
  rank_difference_score.actual <- dirac.rank_difference_score(expression,
                                                              c(1,3,5,7),
                                                              c(2,4,6,8),
                                                              c(2,4,5,6))
  expect_equal(rank_difference_score.actual, rank_difference_score.expected)
})


testthat::test_that("Classification rate works", {
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create list for gene_index
  gene_index <- c(1,4,5,7,10)
  # Create expected classificaiton rate variable
  classification_rate.expected <- 0.875
  # Find the actual classificaiton rate
  classification_rate.actual <- dirac.classification_rate(gene_index,
                                                          expression,
                                                          c(1,2,3,4),
                                                          c(5,6,7,8))
  expect_equal(classification_rate.actual, classification_rate.expected)
})

testthat::test_that("Classification rate comparison works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create list of gene_index vectors
  gene_index_list <- list(A=c(1,4,5), B=c(2,3,6), C=c(7,8,9,10))
  # Create expected list
  classification_rate_comparison.expected <- c(
    A=0.750, B=0.500, C=0.625
  )
  classification_rate_comparison.serial <- dirac.classification_rate.compare(
    expression, c(1,2,3,4), c(5,6,7,8), gene_index_list, parallel=FALSE,
    as.frame=FALSE)
  classification_rate_comparison.parallel <- dirac.classification_rate.compare(
    expression, c(1,2,3,4), c(5,6,7,8), gene_index_list, parallel=TRUE,
    cores=2,as.frame=FALSE)
  expect_equal(classification_rate_comparison.serial,
               classification_rate_comparison.expected)
  expect_equal(classification_rate_comparison.parallel,
               classification_rate_comparison.expected)

})

testthat::test_that("Classification rate best works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create list of gene_index vectors
  gene_index_list <- list(A=c(1,4,5), B=c(2,3,6), C=c(7,8,9,10))
  best_classification_rate.expected <- "A"
  best_classification_rate.actual <- dirac.classification_rate.best(
    expression, c(1,2,3,4), c(5,6,7,8), gene_index_list, parallel=FALSE)
  expect_equal(best_classification_rate.actual,
               best_classification_rate.expected)
})

testthat::test_that("Single compare phenotype works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create the rank matrix
  rank_matrix <- dirac.rank_matrix(expression)
  abs_diff.expected <- 0.005555556
  abs_diff.actual <- dirac.compare_phenotype.single(rank_matrix,
                                                    c(1,2,3,4),
                                                    c(5,6,7,8))
  expect_equal(abs_diff.actual, abs_diff.expected, tolerance=1e-06)
})

testthat::test_that("Shuffled compare phenotype works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  # Create the rank matrix
  rank_matrix <- dirac.rank_matrix(expression)
  abs_diff_comp.expected <- 0.08333333
  abs_diff_comp.actual <- dirac.compare_phenotype.shuffle(
    10, rank_matrix, c(1,2,3,4,5,6,7,8), 4,4, TRUE
  )
  expect_equal(abs_diff_comp.actual, abs_diff_comp.expected, tolerance=1e-06)
})

testthat::test_that("Compare phenotype works",{
  # Create test expression set
  set.seed(42)
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  compare_phenotype_serial.expected <- list(value=0, pval=0.8)
  compare_phenotype_serial.actual <- dirac.compare_phenotype(
    c(1,2,3,4), expression, c(1,2,3,4), c(5,6,7,8), bootstrap_iterations = 20,
    parallel=FALSE, replace=TRUE, seed=42
  )
  expect_equal(compare_phenotype_serial.actual,
               compare_phenotype_serial.expected)
  compare_phenotype_parallel.expected <- list(value=0, pval=0.9)
  compare_phenotype_parallel.actual <- dirac.compare_phenotype(
    c(1,2,3,4), expression, c(1,2,3,4), c(5,6,7,8), bootstrap_iterations = 20,
    parallel=TRUE, cores=2, replace=TRUE, seed=42
  )
  expect_equal(compare_phenotype_parallel.actual,
               compare_phenotype_parallel.expected)
})














