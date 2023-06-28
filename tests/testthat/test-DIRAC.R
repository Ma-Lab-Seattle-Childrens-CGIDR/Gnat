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
  template.actual <- dirac.rank_template(expression)
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
