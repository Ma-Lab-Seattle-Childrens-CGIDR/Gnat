testthat::test_that("Entropy calculation works for vectors", {
  test_vector <- c(1,1,1,1,1)
  testthat::expect_equal(entropy.vector(test_vector), 0)
  test_vector2 <- c(1,1,2)
  testthat::expect_equal(entropy.vector(test_vector2), 0.9182958,
                         tolerance = 1e-4)
  test_vector3 <- c(1,2,3)
  testthat::expect_equal(entropy.vector(test_vector3), 1.584963,
                         tolerance = 1e-4)
})

testthat::test_that("Entropy calculation works for matrices",{
  test_matrix <- matrix(c(1,1,1,1,1,
                          1,2,1,2,2,
                          1,2,3,4,5,
                          6,5,4,3,2), ncol=4)
  expect_length(entropy.matrix(test_matrix, margin=2), 4)
  expect_length(entropy.matrix(test_matrix, margin=1), 5)
  entropy.expected <- c(0,0.9709506, -log2(1/5), -log2(1/5))
  expect_equal(entropy.matrix(test_matrix, margin=2),
               entropy.expected)
})

testthat::test_that("Computing rank matrix works",{
  test_matrix <- matrix(c(0.76,1.98,2.5,4.4,5.1,
                          5.6,1.2,2.2,4.5,3.1,
                          1.05,4.09,2.9,3.2,5.3,
                          2.3,3.2,1.4,4.5,5.4), ncol=4)
  expected_matrix <- matrix(c(1,2,3,4,5,
                              5,1,2,4,3,
                              1,4,2,3,5,
                              2,3,1,4,5), ncol=4)
  rank_matrix.actual <- infer.rank_matrix(test_matrix)
  expect_equal(nrow(rank_matrix.actual), nrow(expected_matrix))
  expect_equal(ncol(rank_matrix.actual), ncol(expected_matrix))
  expect_equal(rank_matrix.actual, expected_matrix)
})

testthat::test_that("Gene entropy works", {
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  gene_entropy.actual <-  infer.gene.entropy(expression, c(1,5,6,7,9),
                                             c(1,2,3,4))
  gene_entropy.expected <- c(1.5, 1.0, 0.8112781, 1.5, 0.8112781)
  expect_length(gene_entropy.actual, 5)
  expect_equal(gene_entropy.actual, gene_entropy.expected, tolerance=1e-6)
})

testthat::test_that("Gene network entropy works", {
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  gene_network_entropy.actual <-
    infer.gene_network.entropy(expression, gene_network = c(1,2,3,4,5,6,7,8),
                               phenotype = c(1,2,3,4,5,6,7,8))
  gene_network_entropy.expected <-
    mean(infer.gene.entropy(expression, gene_network = c(1,2,3,4,5,6,7,8),
                            phenotype = c(1,2,3,4,5,6,7,8)))
  expect_equal(gene_network_entropy.actual, gene_network_entropy.expected)
})

testthat::test_that("Single compare phenotypes works", {
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  abs_diff.actual <- infer.compare_phenotypes.single(
    infer.rank_matrix(expression[c(1,2,3,6,7),]),c(1,2,3,4), c(5,6,7,8))
  abs_diff.expected <- 0.53774438
  expect_equal(abs_diff.actual, abs_diff.expected, tolerance=1e-6)
})

testthat::test_that("Shuffled compare phenotypes works", {
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  abs_diff.actual <-
    infer.compare_phenotypes.shuffle(
      10, infer.rank_matrix(expression[c(1,2,3,4,5,6),]), c(1,2,3,4,5,6,7,8),
      4, 4, replace = TRUE)
  expect_type(abs_diff.actual, "double")
  abs_diff.expected <- 0.135213
  expect_equal(abs_diff.actual, abs_diff.expected, tolerance = 1e-5)
})

testthat::test_that("Compare phenotypes works", {
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  gene_network <- c(1,2,3,4,5,6)
  phenotype1 <- c(1,2,3,4)
  phenotype2 <- c(5,6,7,8)
  bootstrap_iterations <- 1000
  seed <- 42
  os_type <- .Platform$OS.type
  actual.parallel <- infer.compare_phenotypes(
    c(1,2,3,4,5,6), expression, c(1,2,3,4), c(5,6,7,8),
    bootstrap_iterations = 100, parallel=TRUE, cores=2, replace=TRUE, seed=42)
  expect_named(actual.parallel, c("p1.rank_entropy", "p2.rank_entropy",
                                  "absolute_difference", "p.value"))
  expect_gt(actual.parallel$p.value, 0.10)
  expect_equal(abs(actual.parallel$p1.rank_entropy -
                     actual.parallel$p2.rank_entropy),
               actual.parallel$absolute_difference)



  expression2.p1 <- matrix(c(1,2,3,4,5,6,7,8,
                          1,2,3,4,5,6,8,7,
                          1,2,3,4,5,6,7,8,
                          1,2,3,4,5,7,6,8,
                          1,3,2,4,5,6,7,8,
                          1,2,3,4,5,6,7,8,
                          1,2,3,4,6,5,7,8), ncol=7)
  expression2.p2 <- matrix(runif(7*8, min=1, max=500), ncol = 7)
  expression2 <- cbind(expression2.p1, expression2.p2)
  actual2.parallel <- infer.compare_phenotypes(
    c(1,2,3,4,5,6,7,8), expression2, c(1,2,3,4,5,6,7), c(8,9,10,11,12,13,14),
    bootstrap_iterations = 100, parallel=TRUE, cores=2, replace=TRUE, seed=42
  )
  expect_named(actual2.parallel, c("p1.rank_entropy", "p2.rank_entropy",
                                  "absolute_difference", "p.value"))
  expect_lt(actual2.parallel$p.value, 0.05)
  expect_equal(abs(actual2.parallel$p1.rank_entropy -
                     actual2.parallel$p2.rank_entropy),
               actual2.parallel$absolute_difference)

  # Test serial operation
  actual.serial <- infer.compare_phenotypes(
    c(1,2,3,4,5,6), expression, c(1,2,3,4), c(5,6,7,8),
    bootstrap_iterations = 100, parallel=FALSE, cores=2, replace=TRUE, seed=42)
  expect_named(actual.serial, c("p1.rank_entropy", "p2.rank_entropy",
                                   "absolute_difference", "p.value"))
  expect_gt(actual.serial$p.value, 0.10)
  expect_equal(abs(actual.serial$p1.rank_entropy -
                     actual.serial$p2.rank_entropy),
               actual.serial$absolute_difference)
})

testthat::test_that("INFER works",{
  os_type <- .Platform$OS.type
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  gene_network_list <- list(A=c(1,2,3), B=c(4,5,6), C=c(7,8,9,10))
  p1 <- c(1,2,3,4)
  p2 <- c(5,6,7,8)
  res.actual <- INFER(expression, gene_network_list = gene_network_list,
                      phenotype1 = p1, phenotype2 = p2,
                      bootstrap_iterations = 100, parallel=TRUE, cores=2,
                      replace=TRUE, seed=42, as.frame=TRUE)
  expect_equal(nrow(res.actual),3)
  expect_equal(ncol(res.actual), 5)
  expect_equal(colnames(res.actual), c("gene_network", "p1.rank_entropy",
                                       "p2.rank_entropy",
                                         "absolute_difference", "p.value"))
  if(os_type=="windows"){
    res.expected <- data.frame(
      gene_network = c("A","B","C"),
      p1.rank_entropy = c(1.1037594,1.3333333,0.8278195),
      p2.rank_entropy = c(1.103759, 1.103759, 1.452820),
      absolute_difference = c(0.000000, 0.229574, 0.625000),
      p.value = c(0.77,0.64, 0.22)
    )
  } else{
    res.expected <- data.frame(
      gene_network = c("A","B","C"),
      p1.rank_entropy = c(1.1037594,1.3333333,0.8278195),
      p2.rank_entropy = c(1.103759, 1.103759, 1.452820),
      absolute_difference = c(0.000000, 0.229574, 0.625000),
      p.value = c(0.79,0.62, 0.19)
    )
  }
  rownames(res.expected) <- c("A","B","C")
  expect_equal(colnames(res.expected), c("gene_network", "p1.rank_entropy",
                                         "p2.rank_entropy",
                                         "absolute_difference", "p.value"))
  expect_equal(res.actual, res.expected, tolerance = 1e-5)

  # Test serial operation
  res.actual <- INFER(expression, gene_network_list = gene_network_list,
                      phenotype1 = p1, phenotype2 = p2,
                      bootstrap_iterations = 100, parallel=FALSE, cores=1,
                      replace=TRUE, seed=42, as.frame=TRUE)
  expect_equal(nrow(res.actual),3)
  expect_equal(ncol(res.actual), 5)
  expect_equal(colnames(res.actual), c("gene_network", "p1.rank_entropy",
                                         "p2.rank_entropy",
                                         "absolute_difference", "p.value"))
  expect_equal(
    abs(res.actual$p1.rank_entropy-res.actual$p2.rank_entropy),
    res.actual$absolute_difference
  )
  # Difference due to different RNG between serial and parallel
  if(os_type=="windows"){
    res.expected <- data.frame(
      gene_network = c("A","B","C"),
      p1.rank_entropy = c(1.1037594,1.3333333,0.8278195),
      p2.rank_entropy = c(1.103759, 1.103759, 1.452820),
      absolute_difference = c(0.000000, 0.229574, 0.625000),
      p.value = c(0.71,0.54, 0.15)
    )
  } else {
    res.expected <- data.frame(
      gene_network = c("A","B","C"),
      p1.rank_entropy = c(1.1037594,1.3333333,0.8278195),
      p2.rank_entropy = c(1.103759, 1.103759, 1.452820),
      absolute_difference = c(0.000000, 0.229574, 0.625000),
      p.value = c(0.71,0.54, 0.15)
    )
  }
  rownames(res.expected) <- c("A","B","C")
  expect_equal(res.actual, res.expected, tolerance = 1e-5)
})








