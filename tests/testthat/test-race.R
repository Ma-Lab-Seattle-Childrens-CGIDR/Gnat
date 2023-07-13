test_that("Mean rank correlation works",{
  test_matrix <- matrix(1:16, ncol=4)
  testthat::expect_equal(race.mean_rank_corr(test_matrix), 1)
  test_matrix2 <- matrix(c(1,2,3,4,1,2,4,3,1,4,2,3,4,3,2,1), ncol=4)
  mean_rank_correlation <- race.mean_rank_corr(test_matrix2)
  expect_lte(mean_rank_correlation, 1)
  expect_gte(mean_rank_correlation, -1)
  expect_equal(mean_rank_correlation, -0.166666667, tolerance=1e-5)
})

test_that("Compare single phenotypes works",{
  set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
  expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
  result <- race.compare_phenotypes.single(expression[c(1,2,3,4),], c(1,2,3,4),
                                           c(5,6,7,8))
  expect_lte(result, 2)
  expect_gte(result, -2)
  expect_equal(result, 0.277777778, tolerance=1e-5)
})

test_that("Compare shuffled phenotypes works", {
   set.seed(42, kind = "Mersenne-Twister", normal.kind = "Inversion")
   expression <- matrix(floor(runif(10*8, min=1, max=500000)), ncol=8, nrow=10)
   result <- race.compare_phenotypes.shuffle(10, expression[c(1,4,5,6,8),],
                                             c(1,2,3,4,5,6,7,8), 4,4,
                                             replace=TRUE)
   expect_lte(result,2)
   expect_gte(result,-2)
   expect_equal(result, 0.16666667, tolerance = 1e-5)
})

test_that("Compare phenotypes works",{
  expression.p1 <- matrix(c(1,2,3,4,5,6,7,8,
                             1,2,3,4,5,6,8,7,
                             1,2,3,4,5,6,7,8,
                             1,2,3,4,5,7,6,8,
                             1,3,2,4,5,6,7,8,
                             1,2,3,4,5,6,7,8,
                             1,2,3,4,6,5,7,8), ncol=7)
  expression.p2 <- matrix(runif(7*8, min=1, max=500), ncol = 7)
  expression <- cbind(expression.p1, expression.p2)
  result.serial <- race.compare_phenotypes(c(1,3,4,6,7), expression,
                                    phenotype1 = 1:7, phenotype2=8:14,
                                    bootstrap_iterations = 100,
                                    parallel=FALSE, cores=1, replace=TRUE,
                                    seed=42)
  result.parallel <- race.compare_phenotypes(c(1,3,4,6,7), expression,
                                           phenotype1 = 1:7, phenotype2=8:14,
                                           bootstrap_iterations = 100,
                                           parallel=TRUE, cores=2,
                                           replace=TRUE, seed=42)
  expect_named(result.serial, c("p1.mean_rank_correlation",
                                "p2.mean_rank_correlation",
                                "absolute_difference",
                                "p.value"))
  expect_named(result.parallel, c("p1.mean_rank_correlation",
                                "p2.mean_rank_correlation",
                                "absolute_difference",
                                "p.value"))
  expect_lte(result.serial$p.value, 0.05)
  expect_lte(result.serial$p.value, 0.05)
  random_expression <- matrix(rnorm(8*10), ncol=8)
  result.serial <- race.compare_phenotypes(c(1,3,5,6,8),
                                           random_expression, c(1,2,3,4),
                                           c(5,6,7,8),
                                           bootstrap_iterations = 100,
                                           parallel=FALSE, cores=1,
                                           replace=TRUE, seed=42)
  result.parallel <- race.compare_phenotypes(c(1,3,5,6,8),
                                           random_expression, c(1,2,3,4),
                                           c(5,6,7,8),
                                           bootstrap_iterations = 100,
                                           parallel=TRUE, cores=2,
                                           replace=TRUE, seed=42)
  expect_gt(result.serial$p.value, 0.05)
  expect_gt(result.parallel$p.value, 0.05)
})

test_that("RACE works",{
  expression.p1 <- matrix(c(1,2,3,4,5,6,7,8,
                            1,2,3,4,5,6,8,7,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,5,7,6,8,
                            1,3,2,4,5,6,7,8,
                            1,2,3,4,5,6,7,8,
                            1,2,3,4,6,5,7,8), ncol=7)
  expression.p2 <- matrix(runif(7*8, min=1, max=500), ncol = 7)
  expression <- cbind(expression.p1, expression.p2)
  gene_network_list <- list(A=c(1,4,7,8), B=c(2,3,6), C=c(2,4,5))
  result.serial <- RACE(expression, gene_network_list,
                        phenotype1=1:7, phenotype2 = 8:14,
                        bootstrap_iterations = 100,
                        parallel = FALSE,cores=1, replace=TRUE, seed=42,
                        as.frame=TRUE)
  result.parallel <- RACE(expression, gene_network_list,
                          1:7, 8:14, bootstrap_iterations = 200,
                          parallel=TRUE, core=2, replace=TRUE, seed=42,
                          as.frame=TRUE)
  expect_equal(colnames(result.serial), c("gene_network",
                                          "p1.mean_rank_correlation",
                                          "p2.mean_rank_correlation",
                                          "absolute_difference",
                                          "p.value"))
  expect_equal(colnames(result.parallel), c("gene_network",
                                          "p1.mean_rank_correlation",
                                          "p2.mean_rank_correlation",
                                          "absolute_difference",
                                          "p.value"))
  expect_equal(result.serial$gene_network, c("A","B","C"))
  expect_equal(result.parallel$gene_network, c("A","B","C"))
  expect_equal(result.serial$p1.mean_rank_correlation,
               result.parallel$p1.mean_rank_correlation)
  expect_equal(result.serial$p2.mean_rank_correlation,
               result.parallel$p2.mean_rank_correlation)
  expect_true(all(result.serial$p.value<=0.05))
  expect_true(all(result.parallel$p.value<=0.05))
})





