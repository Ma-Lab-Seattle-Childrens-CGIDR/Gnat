test_that("Rank matrix works",{
  test_matrix <- matrix(1:16, ncol=4)
  rank_matrix.actual <- crane.rank_matrix(test_matrix)
  expect_equal(nrow(rank_matrix.actual),4)
  expect_equal(ncol(rank_matrix.actual),4)
  rank_matrix.expected <- matrix(rep(c(1,2,3,4), 4), ncol=4)
  expect_equal(rank_matrix.actual, rank_matrix.expected)
})

test_that("Mean centroid distance works", {
  test_matrix <- matrix(1:16, ncol=4)
  expect_equal(crane.mean_centroid_distance(crane.rank_matrix(test_matrix)), 0)
  set.seed(42, kind="Mersenne-Twister", normal.kind="Inversion")
  test_matrix <- matrix(rnorm(16), ncol=4)
  result <- crane.mean_centroid_distance(
    crane.rank_matrix(test_matrix)
  )
  expect_equal(result, 1.620591, tolerance=1e-5)
})

test_that("Single compare phenotypes works",{
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  rank_matrix <- crane.rank_matrix(expression)
  result <- crane.compare_phenotypes.single(rank_matrix = rank_matrix,
                                            c(1,2,3,4), c(5,6,7,8))
  expect_type(result, "double")
  expect_equal(result, 0.1458526, tolerance=1e-5)
  result2 <- crane.compare_phenotypes.single(rank_matrix,
                                             c(1,2,3,4),
                                             c(1,2,3,4))
  expect_equal(result2, 0)
})

test_that("Shuffle compare phenotypes works",{
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  rank_matrix <- crane.rank_matrix(expression)
  result <- crane.compare_phenotypes.shuffle(10, rank_matrix,
                                             1:8,4,4,TRUE)
  expect_type(result, "double")
  expect_equal(result, 1.439112, tolerance=1e-5)
})

test_that("Compare Phenotypes works", {
  expression.p1 <- matrix(c(1,2,3,4,5,6,7,8,
                            1,2,3,4,5,6,8,7,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,5,7,6,8,
                            1,3,2,4,5,6,7,8,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,6,5,7,8), ncol=7)
  expression.p2 <- matrix(runif(7*8, min=1, max=500), ncol = 7)
  expression <- cbind(expression.p1, expression.p2)
  result.serial <- crane.compare_phenotypes(c(1,3,4,6,7), expression,
                                           phenotype1 = 1:7, phenotype2=8:14,
                                           bootstrap_iterations = 100,
                                           parallel=FALSE, cores=1, replace=TRUE,
                                           seed=42)
  result.parallel <- crane.compare_phenotypes(c(1,3,4,6,7), expression,
                                             phenotype1 = 1:7, phenotype2=8:14,
                                             bootstrap_iterations = 100,
                                             parallel=TRUE, cores=2,
                                             replace=TRUE, seed=42)
  expect_named(result.serial, c("p1.mean_centroid_distance",
                                "p2.mean_centroid_distance",
                                "absolute_difference",
                                "p.value"))
  expect_named(result.parallel, c("p1.mean_centroid_distance",
                                "p2.mean_centroid_distance",
                                "absolute_difference",
                                "p.value"))
  expect_lte(result.serial$p.value, 0.1)
  expect_lte(result.parallel$p.value, 0.1)
  random_expression <- matrix(rnorm(8*10), ncol=8)
  result.random <- crane.compare_phenotypes(c(1,2,3,4,10), random_expression,
                                            phenotype1=1:4, phenotype2=5:8,
                                            bootstrap_iterations = 100,
                                            parallel=FALSE, cores=1,
                                            replace=TRUE, seed=42)
  expect_gte(result.random$p.value, 0.1)
})

test_that("CRANE works",{
  expression.p1 <- matrix(c(1,2,3,4,5,6,7,8,
                            1,2,3,4,5,6,8,7,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,5,7,6,8,
                            1,3,2,4,5,6,7,8,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,6,5,7,8), ncol=7)
  expression.p2 <- matrix(runif(7*8, min=1, max=500), ncol = 7)
  expression <- cbind(expression.p1, expression.p2)
  gene_network_list <- list(A=c(1,2,3), B=c(4,5,6), C=c(1,5,8))
  p1=1:7
  p2=8:14
  result <- CRANE(expression, gene_network_list, p1, p2,
                  bootstrap_iterations=100, parallel=TRUE, cores=2,
                  seed=42, as.frame = TRUE)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
  expect_equal(colnames(result), c("gene_network",
                                    "p1.mean_centroid_distance",
                                    "p2.mean_centroid_distance",
                                    "absolute_difference",
                                    "p.value"))
  expect_equal(result$gene_network, c("A","B","C"))
  expect_true(all(result$p.value<=0.1))
})
